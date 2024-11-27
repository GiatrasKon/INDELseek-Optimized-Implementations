#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <regex>
#include <limits>
#include <iomanip>
#include <chrono>
#include <unistd.h>
#include <sys/resource.h>

using StringVector = std::vector<std::string>;
using chr_t = uint32_t;
using position_t = uint32_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Helper functions for RAM logging during execution
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double getRAMUsage() {
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) != 0) {
    return -1.0; // Error occurred
  }

  return static_cast<double>(usage.ru_maxrss); // Convert from KB to MB
}

void logRamUsage(std::ofstream& logfile) {
//  std::ofstream logfile("ram_usage.log", std::ios::app);
  if (not logfile.is_open()) {
    std::cerr << "Failed to open log file." << std::endl;
    exit(1);
  }

  // Get current time
  time_t currentTime = time(nullptr);
  tm* localTime = localtime(&currentTime);
  char timeString[80];
  strftime(timeString, 80, "%Y-%m-%d %H:%M:%S", localTime);

  // Log the RAM usage and current time to the file
  logfile << "[" << timeString << "] RAM Usage: " << getRAMUsage() << " kB" << std::endl;
  logfile.flush();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Class that represents the INDELSeek algorithm, with read BAM alignments capabilities,
/// complex indel detection and write methods to VCF.
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class IndelSeek {
public:
  IndelSeek() : _ramLogFile("ram_usage_" + std::to_string(getpid()) + ".log", std::ios::app) {}
  ~IndelSeek() { _ramLogFile.close(); }
  void parseArguments(int argC, char* argV[]);
  void printOptions();
  void processFile();
  void writeMutationsInVCF();

private:
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Struct representing a CIGAR Operation, like 10M (10 base matches)
/// It's parsed so we have easy access to length and operation character,
/// along with compare methods for use in hash maps.
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  struct CigarOp {
    unsigned long _length = 0;
    char     _operation;

    bool operator==(CigarOp const& obj) const {
      return _length == obj._length and _operation == obj._operation;
    }

    bool operator!=(CigarOp const& obj) const {
      return not (*this == obj);
    }
  };

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Struct representing a CIGAR sequence, like 10M5D3M20X
/// It's parsed so we have easy access to each operation via indexing
/// and a function that calculates the total length.
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  struct CigarSequence {
    std::vector<CigarOp> _sequence;
    size_t totalLength() const;
    void   print() const;
  };

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Basically a representation of a possible mutation, start and end coordinates,
/// reference and variant sequences and the average of all position qualities.
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  struct AlignmentCluster {
    uint32_t _start = 0;
    uint32_t _end = 0;
    uint32_t _minPos = 0;
    uint32_t _maxPos = 0;
    std::string _refSeq;
    std::string _readSeq;
    float _meanQual = 0;
  };

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Those stats are to measure the appearance of each mutation
/// in the whole BAM file (positive or negative strand).
/// Used to calculate the allelic frequency.
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  struct MutationStats {
    size_t _positiveStrand = 0, _negativeStrand = 0;
    std::vector<float> _positiveStrandQual, _negativeStrandQual;
  };

  CigarSequence parseCigar(const std::string& cigar);
  std::vector<AlignmentCluster> reconstructAlignment(const CigarSequence& cigarSequence, const std::string& readSeq, const std::string& readQual, const std::string& refSeq, const std::string& startPosStr);
  MutationStats* findMutationInCache(const std::string& chromosome, const position_t start, const std::string& key);
  std::string faidx(const std::string &genomicRegion);
  static std::string regionRepresentation(const std::string& chromosome, const std::string& start, const CigarSequence& cigar);
  static void updateStats(MutationStats& stats, const std::string& direction, float meanQual);
  double depth(const std::string& genomicRegion);
  void printVCFHeader();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Thresholds and flags for the various filters that used to filter false positive detections
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  uint32_t _phredOffset = 33,
           _qualityThreshold = 20,
           _minDepth = 50,
           _maxSamtoolsDepth = 500000,
           _maxDistance = 5;

  double   _minAF = 0;
  bool     _skipLowQual = false,
           _skipLowDepth = false,
           _skipLowAF = false;

  /// the various filenames used
  std::string _refGenomeFile = "ucsc.hg19.fasta", _inputFileName, _outputFileName, _depthBam;
  std::ofstream _ramLogFile;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// The main storing structure for the detected mutations.
/// It's a threeway hashMap storing each mutation per position found:
/// chromosome -> position -> mutationKey -> mutation statistics
///
/// Additional advantage is that map (in C++) is de facto sorted. So when we write into VCF we don't sort explicitly.
///
/// The mutation key is a unique string created for each mutation: chromosome|position|refseq|readseq
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::map<chr_t, std::map<position_t, std::map<std::string, MutationStats>>> _mutations;
  size_t totalMutations();
};

/// Calculates the total number of mutation that we have found
size_t IndelSeek::totalMutations()
{
  size_t total = 0;
  for (const auto& [chrNum, positions] : _mutations) {
    for (const auto &[startNum, keys]: positions) {
      total += keys.size();
    }
  }

  return total;
}

/// Helper function for concatinating strings with specific delimiter
template<typename StringType, typename Container>
std::string join(const StringType& separator, const Container& c)
{
  if (c.size() == 0)
    return {};
  std::stringstream s;
  auto i = c.begin();
  s << *i;
  for (++i; i != c.end(); ++i)
    s << separator << *i;
  return s.str();
}

/// Checks if a string starts with a specific substring
bool startsWith(const std::string& text, const std::string& prefix)
{
  return text.rfind(prefix, 0) == 0;
}

/// Total length of a CIGAR Sequence
size_t IndelSeek::CigarSequence::totalLength() const
{
  size_t result = 0;
  for (const auto& cig : _sequence)
    result += cig._length;

  return result;
}

void IndelSeek::CigarSequence::print() const
{
  for (const auto& cig : _sequence)
    std::cout << cig._length << cig._operation << (cig != _sequence.back() ? "-" : "");
  std::cout << std::endl;
}

/// Concatinates string representing a region for search in Samtools
std::string IndelSeek::regionRepresentation(const std::string &chromosome,
                                            const std::string &start,
                                            const CigarSequence &cigar)
{
  std::string region = chromosome + ":" + start + "-" + std::to_string(std::stoul(start) + cigar.totalLength() - 1);
  if (not startsWith(chromosome,"chr"))
    region = "chr" + region;

  return region;
}

/// Parses a CIGAR string into each operation & length, returning a dynamic array (vector)
IndelSeek::CigarSequence IndelSeek::parseCigar(const std::string &cigar)
{
  CigarSequence sequence;

  std::string numStr;
  for (const auto c : cigar) {
    if (std::isdigit(c))
      numStr += c;
    else {
      if (c == 'N' or c == 'P') {
        std::cout << "ERROR: Unexpected CIGAR operation: " << c << std::endl;
        exit(1);
      }

      if (c != 'H') {
        try {
          CigarOp cigarOp = { std::stoul(numStr), c };
          sequence._sequence.push_back(cigarOp);
        }
        catch (const std::out_of_range& e) {
          std::cout << "ERROR: Can't convert number " << numStr << " to unsigned integer." << std::endl;
          throw;
        }
      }

      numStr = "";
    }
  }

  return sequence;
}

/// Splits a string into its parts based on the provided delimiter
StringVector split(const std::string& text, const char delimiter = '\t')
{
  StringVector result;
  std::string token;
  std::stringstream stream(text);
  while (getline(stream, token, delimiter))
    result.push_back(token);

  return result;
}

/// Method that parses the command line arguments passed to the program
void IndelSeek::parseArguments(int argC, char* argV[])
{
  const StringVector arguments(argV + 1, argV + argC);

  for (auto i = 0; i < arguments.size(); i++) {
    const auto& argument = arguments[i];
    if (argument.front() == '-') {
      if (argument == "--phred_offset")
        _phredOffset = std::stoul(arguments[++i]);
      if (argument == "--qual_threshold")
        _qualityThreshold = std::stoul(arguments[++i]);
      if (argument == "--min_depth")
        _minDepth = std::stoul(arguments[++i]);
      if (argument == "--max_depth")
        _maxSamtoolsDepth = std::stoul(arguments[++i]);
      if (argument == "--max_distance")
        _maxDistance = std::stoul(arguments[++i]);

      if (argument == "--reference")
        _refGenomeFile = argument[++i];

      if (argument == "--input_sam")
        _inputFileName = arguments[++i];

      if (argument == "--output_vcf")
        _outputFileName = arguments[++i];

      if (argument == "--depth_bam")
        _depthBam = arguments[++i];

      if (argument == "--min_af")
        _minAF = std::stod(arguments[++i]);

      if (not _skipLowQual)
        _skipLowQual = (argument == "--skip_low_qual");

      if (not _skipLowDepth)
        _skipLowDepth = (argument == "--skip_low_depth");

      if (not _skipLowAF)
        _skipLowAF = (argument == "--skip_low_af");
    }
  }

  if (_outputFileName.empty()) {
    std::cout << "ERROR: --output_vcf must be set." << std::endl;
    exit(1);
  }

  if (_inputFileName.empty())
    std::cout << "Input file was empty so reading from standard input." << std::endl;
}

std::string boolToYesNo(bool value)
{
  return value ? "YES" : "NO";
}

/// A method to output all the thresholds and options into standard output and log file.
void IndelSeek::printOptions()
{
  std::stringstream ss;
  ss << "\tInput BAM file: " << _inputFileName << std::endl
     << "\tOutput VCF: " << _outputFileName << std::endl
     << "\tReference genome: " << _refGenomeFile << std::endl
     << "\tPhred Offset=" << _phredOffset << std::endl
     << "\tQuality Threshold=" << _qualityThreshold << std::endl
     << "\tMin Depth=" << _minDepth << std::endl
     << "\tMax Depth=" << _maxSamtoolsDepth << std::endl
     << "\tMax Distance=" << _maxDistance << std::endl
     << "\tDepth BAM (used for calculating AF)=" << _depthBam << std::endl
     << "\tMin AF=" << _minAF << std::endl
     << "\tSkip Low Qual=" << boolToYesNo(_skipLowQual) << std::endl
     << "\tSkip Low Depth=" << boolToYesNo(_skipLowDepth) << std::endl
     << "\tSkip Low AF=" << boolToYesNo(_skipLowAF) << std::endl;

  std::cout << ss.str() << std::endl;
  _ramLogFile << ss.str() << std::endl; //so we know from which run the log is
}

/// Method that executes a command directly to the system and returns the output of that command
std::string exec(const std::string& cmd) {
  std::array<char, 128> buffer{};
  std::string result;
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
  }

  return result;
}

/// Function for calculating the mean value of a sequential container of values
template<typename T>
double mean(const std::vector<T>& container)
{
  T sum = 0;
  for (const auto entry : container)
    sum += entry;

  return (double) sum / (double) container.size();
}

/// Method that calls samtools to get the reference sequence of a given region
static size_t regionCacheHits = 0;
std::string IndelSeek::faidx(const std::string &genomicRegion)
{
  /// The cache (hash map) is used to store the refseq of regions already queried, so we don't call samtools again
  static std::map<std::string, std::string> cache;
  const auto found = cache.find(genomicRegion);
  if (found != cache.end()) {
    regionCacheHits++;
    return found->second;
  }


  const auto command = "samtools faidx " + _refGenomeFile + " " + genomicRegion;
  const auto& output = exec(command);

  /// Parsing the output of samtools faidx command
  const auto entries = split(output, '\n');
  if (entries[0] != (">" + genomicRegion)) {
    std::cout << "ERROR: Faidx region: " << entries[0] << " doesn't match requested region: " << genomicRegion << std::endl;
    exit(1);
  }

  std::string sequence;
  for (const auto& entry : entries) {
    if (entry == entries.front())
      continue;

    /// capitalise all the letters just to be sure that a == A, t == T etc.
    for (const auto& c : entry)
      sequence += std::toupper(c);
  }

  /// store in cache for future use
  cache[genomicRegion] = sequence;
  return sequence;
}

/// Method to query the depth of a given region from a BAM file (depth == times of a specific sequence existing in BAM)
static size_t depthCacheHits = 0;
double IndelSeek::depth(const std::string& genomicRegion)
{
  /// The cache (hash map) is used to store the depth of regions already queried, so we don't call samtools again
  static std::map<std::string, double> cache;
  const auto found = cache.find(genomicRegion);
  if (found != cache.end()) {
    depthCacheHits++;
    return found->second;
  }

  const auto command = "samtools depth -d " + std::to_string(_maxSamtoolsDepth) +
                                  " -r " + genomicRegion +
                                  " -q " + std::to_string(_qualityThreshold) +
                                  " " + _depthBam;

  const auto& output = exec(command);
  std::vector<uint32_t> depths;
  /// Parsing the depth for each position in the region
  for (const auto& line : split(output, '\n')) {
    const auto lineParts = split(line);
    if (lineParts.size() != 3) {
      std::cout << "ERROR: Invalid depth=" << line << " received for region " << genomicRegion << std::endl;
      exit(1);
    }
    depths.emplace_back(std::stoul(lineParts[2]));
  }

  /// Mean depth of all the positions, represents the mean depth of the region
  const auto meanDepths = mean(depths);

  /// Store in cache
  cache[genomicRegion] = meanDepths;

  return meanDepths;
}

/// Find the minimum value in a vector (dynamic array) of strings (string representation of numbers coming from reading the BAM)
uint32_t minFromStringVector(const StringVector& vectorOfStrings, int start, int end)
{
  uint32_t min = std::numeric_limits<uint32_t>::max();
  for (int i = start; i <= end; i++) {
    const auto& entry = vectorOfStrings[i];
    if (entry == "*" or entry.empty())
      continue;

    const auto number = (uint32_t)std::stoul(entry);
    min = std::min(min, number);
  }
  return min;
}

/// Find the maximum value in a vector (dynamic array) of strings (string representation of numbers coming from reading the BAM)
uint32_t maxFromStringVector(const StringVector& vectorOfStrings, int start, int end)
{
  int max = -1;
  for (int i = start; i <= end; i++) {
    const auto& entry = vectorOfStrings[i];
    if (entry == "*" or entry.empty())
      continue;

    const auto number = std::stoi(entry);
    max = std::max(max, number);
  }
  return (uint32_t)max;
}

/// Removes * from sequence
std::string cleanSeq(const std::vector<char>& sequence, int start, int end)
{
  std::string result;
  for (int i = start; i <= end; i++) {
    const auto& c = sequence[i];
    if (c == '*')
      continue;

    result += c;
  }
  return result;
}

/// Calculates the mean quality from the quality sequence. Since we need one letter representation for quality
/// (so there is one to one matching with read sequence of nucleotides),
/// ASCII representation is used, for example a quality of A means 65 and a quality of [ means 91.
/// https://web.alfredstate.edu/faculty/weimandn/miscellaneous/ascii/ascii_index.html
float meanQuality(const std::vector<char>& qualseq, int start, int end)
{
  size_t sum = 0;
  for (int i = start; i <= end; i++) {
    const auto& c = qualseq[i];
    if (c == ' ')
      continue;

    sum += (size_t)c; /// here is the conversation to decimal number from ASCII
  }

  return (float)sum/(float)(end-start+1);
}

/// The VCF header for the final output file
void IndelSeek::printVCFHeader()
{
  std::fstream vcf;
  vcf.open(_outputFileName, std::ios::out);

  if (vcf.is_open()) {
    vcf << "##fileformat=VCFv4.1\n";
    vcf << "##source=INDELseek\n";
    vcf << "##reference=file://" << _refGenomeFile << std::endl;
    vcf << "##INFO=<ID=DP2,Number=2,Type=Integer,Description=\"# alt-foward and alt-reverse reads\">\n";
    vcf << "##INFO=<ID=QS2,Number=2,Type=Float,Description=\"Mean quality scores of alt-foward and alt-reverse bases\">\n";
    if (not _depthBam.empty()) {
      vcf << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n";
      vcf << "##INFO=<ID=RDP,Number=1,Type=Float,Description=\"Mean read depth of REF positions\">\n";
      vcf << "##FILTER=<ID=LowAF,Description=\"AF below " << _minAF << "\">\n";
    }
    vcf << "##FILTER=<ID=LowDepth,Description=\"ALT depth below " << _minDepth << "\">\n";
    vcf << "##FILTER=<ID=LowQual,Description=\"Mean quality scores below " << _qualityThreshold << " or ALT contains N\">\n";
    vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"<< std::endl;
  }

  vcf.close();
}

/// check if a substring exists in a given string
bool contains(const std::string& searchText, const std::string& toFind)
{
  return searchText.find(toFind) != std::string::npos;
}

/// The final method that writes the found mutations into the VCF file
void IndelSeek::writeMutationsInVCF()
{
  printVCFHeader();

  std::fstream vcf;
  vcf.open(_outputFileName, std::ios::out | std::ios::app);

  for (const auto& [chrNum, positions] : _mutations) {
    for (const auto& [startNum, keys] : positions) {
      for (const auto& [key, mutationStats] : keys) {

        /// A sanity check that the number of qualities for each mutation matches the numbers actually found in BAM
        if (mutationStats._positiveStrandQual.size() != mutationStats._positiveStrand or
            mutationStats._negativeStrandQual.size() != mutationStats._negativeStrand) {
          std::cout << "ERROR: Invalid mutation stats for key=" << key << std::endl;
          std::cout << "+Qualities=" << mutationStats._positiveStrandQual.size() << ", count=" << mutationStats._positiveStrand << std::endl;
          std::cout << "-Qualities=" << mutationStats._negativeStrandQual.size() << ", count=" << mutationStats._negativeStrand << std::endl;
          exit(1);
        }

        /// Writing the various columns
        const auto mutationParts= split(key, '|');
        size_t index = 0;
        const auto& chromosome = mutationParts[index++];
        const auto& start = mutationParts[index++];
        const auto& end = mutationParts[index++];
        const auto& ref = mutationParts[index++];
        const auto& alt = mutationParts[index++];

        /// Calculate various qualities per strand and total
        const auto qualForward = mean(mutationStats._positiveStrandQual);
        const auto qualReverse = mean(mutationStats._negativeStrandQual);
        const auto combinedQual = (qualForward + qualReverse) / 2;
        const auto combinedDepth = mutationStats._positiveStrand + mutationStats._negativeStrand;

        StringVector filterTags;

        /// Filter for minimum quality
        if (combinedQual < _qualityThreshold or contains(alt, "N")) {
          if (_skipLowQual)
            continue;

          filterTags.emplace_back("LowQual");
        }

        /// Filter for minimum depth
        if (combinedDepth < _minDepth) {
          if (_skipLowDepth)
            continue;

          filterTags.emplace_back("LowDepth");
        }

        std::string info = "DP2=" + std::to_string(mutationStats._positiveStrand) + "," + std::to_string(mutationStats._negativeStrand)
                           + ";QS2=" + std::to_string(qualForward) + "," + std::to_string(qualReverse);

        if (not _depthBam.empty()) {
          /// Query the BAM file for the depth of the region of the mutation
          const auto region = chromosome + ":" + start + "-" + end;
          const auto posDepth = depth(region);
          if (posDepth == 0) {
            if (_skipLowAF)
              continue;

            filterTags.emplace_back("LowAF");
          }
          else {
            /// Calculate allele frequency if depth is not 0
            const auto af = (double)combinedDepth/posDepth;
            if (af < _minAF) {
              if (_skipLowAF)
                continue;

              filterTags.emplace_back("LowAF");
            }
            info += ";AF=" + std::to_string(af) + ";RDP=" + std::to_string(posDepth);
          }
        }

        /// Finaly join all the columns and write the row to VCF output
        vcf << join("\t", StringVector {chromosome, start, ".", ref, alt,
                             std::to_string(combinedQual), (filterTags.empty() ? "PASS" : join(";", filterTags)),
                             info}) << std::endl;
      }
    }
  }
}


/// The main function that detects the possible mutations from a given sequence
std::vector<IndelSeek::AlignmentCluster> IndelSeek::reconstructAlignment(const IndelSeek::CigarSequence &cigarSequence,
                                                                         const std::string &readSeq,
                                                                         const std::string &readQual,
                                                                         const std::string &refSeq,
                                                                         const std::string &startPosStr)
{
  std::vector<AlignmentCluster> output;

  std::vector<char> cigarOp, cigarReadSeq, cigarReadQual, cigarRefSeq;
  StringVector cigarRefPos;

  /// turn the various intermediate storing vectors to 1-based, based on the algorithm
  cigarOp.resize(cigarSequence.totalLength() + 1);
  cigarReadSeq.resize(cigarSequence.totalLength() + 1);
  cigarReadQual.resize(cigarSequence.totalLength() + 1);
  cigarRefPos.resize(cigarSequence.totalLength() + 1);
  cigarRefSeq.resize(cigarSequence.totalLength() + 1);

  size_t cigarPosOffset = 1, readbasePosOffset = 0, refPosOffset = std::stoul(startPosStr), startPos = std::stoul(startPosStr);
  for (const auto& cigar : cigarSequence._sequence) {
    const auto& opLength = cigar._length;
    const auto& op = cigar._operation;

    for (auto i = cigarPosOffset; i < (cigarPosOffset + opLength - 1); i++) {
      if (i >= cigarOp.size() or i >= cigarReadSeq.size() or i>= cigarReadQual.size()) {
        std::cout << "ERROR out of bounds allocation with i=" << i << std::endl;
        std::cout << "Cigar=" << opLength << op << ", readSeq=" << readSeq << ", readQ=" << readQual << ", start=" << startPosStr << std::endl;
        exit(1);
      }

      /// Based on the CIGAR sequence assign the accompanying sequence or metric
      cigarOp[i] = op;
      cigarReadSeq[i] = (op == 'D' ? '*' : readSeq[i + readbasePosOffset]);
      cigarReadQual[i] = (op == 'D' ? '*' : readQual[i + readbasePosOffset]);
      cigarRefPos[i] = (op == 'I' or op == 'S' ? "*" : std::to_string(refPosOffset + i - cigarPosOffset));
      cigarRefSeq[i] = (op == 'I' or op == 'S' ? '*' : refSeq[refPosOffset + i - cigarPosOffset - startPos + 1]);
    }

    if (op == 'D')
      readbasePosOffset -= opLength;
    if (op != 'I' and op != 'S')
      refPosOffset += opLength;
    cigarPosOffset += opLength;
  }

  int start = -1;
  int end = -1 * _maxDistance;
  uint32_t minWindowSize = 2;

  for (size_t i = 1; i <= cigarSequence.totalLength(); i++) {
    if (not (cigarOp[i] == 'M' and (cigarReadSeq[i] != cigarRefSeq[i]) or
        cigarOp[i] == 'I' or
        cigarOp[i] == 'D'))
      continue;

    if ((i - end) > _maxDistance) {
      if (start >= 0) {
        if ((end - start + 1) >= minWindowSize) {
          AlignmentCluster cluster = { (uint32_t)start,
                                       (uint32_t)end,
                                       minFromStringVector(cigarRefPos, start, end),
                                       maxFromStringVector(cigarRefPos, start, end),
                                       cleanSeq(cigarRefSeq, start, end),
                                       cleanSeq(cigarReadSeq, start, end),
                                       meanQuality(cigarReadQual, start, end)
                                       };
          output.emplace_back(cluster);
        }
      }

      start = i;
      end = i;
    }
    else {
      if (i > end)
        end = i;
    }
  }

  if (start >= 0) {
    if ((end - start + 1) >= minWindowSize) {
      AlignmentCluster cluster = { (uint32_t)start,
                                   (uint32_t)end,
                                   minFromStringVector(cigarRefPos, start, end),
                                   maxFromStringVector(cigarRefPos, start, end),
                                   cleanSeq(cigarRefSeq, start, end),
                                   cleanSeq(cigarReadSeq, start, end),
                                   meanQuality(cigarReadQual, start, end)
      };
      output.emplace_back(cluster);
    }
  }

  return output;
}

/// convert chromosome name to number
chr_t convertChromosomeToNum(const std::string& chromosome)
{
  if (startsWith(chromosome, "chr"))
    return std::stoul(chromosome.substr(3, chromosome.length()-3));

  return std::stoul(chromosome);
}

/// Method to check if a mutation was already found and return the stats of that.
IndelSeek::MutationStats* IndelSeek::findMutationInCache(const std::string &chromosome,
                                                        const position_t start,
                                                        const std::string &key)
{
  const chr_t chromosomeNum = convertChromosomeToNum(chromosome);

  /// match the chromosome first
  const auto findChr = _mutations.find(chromosomeNum);
  if (findChr == _mutations.end())
    return nullptr;

  /// match the starting position
  const auto& positions = findChr->second;
  const auto findPos = positions.find(start);
  if (findPos == positions.end())
    return nullptr;

  /// match with the unique mutation key
  const auto& keys = findPos->second;
  const auto findKey = keys.find(key);
  if (findKey == keys.end())
    return nullptr;

  /// returns a pointer to already existed stats if found
  return const_cast<MutationStats *>(&findKey->second);
}

/// Increase the stats according to the strand
void IndelSeek::updateStats(IndelSeek::MutationStats& stats, const std::string& direction, float meanQual)
{
  if (direction == "+") {
    stats._positiveStrand++;
    stats._positiveStrandQual.push_back(meanQual);
  }
  else {
    stats._negativeStrand++;
    stats._negativeStrandQual.push_back(meanQual);
  }
}

/// the main function that reads and process the BAM file
void IndelSeek::processFile()
{
  auto startTime = std::chrono::steady_clock::now();

  /// check if an input file was provided otherwise read from standard input (useful for piping)
  bool readFromInputSAM = false;
  std::fstream sam;
  if (not _inputFileName.empty()) {
    sam.open(_inputFileName, std::ios::in);
    if (not sam.is_open()) {
      std::cout << "ERROR: Can't open file: " << _inputFileName << std::endl;
      exit(1);
    }
    readFromInputSAM = true;
  }

  size_t lineCount = 0, validLines = 0;

  std::string line;
  while (getline(readFromInputSAM ? sam : std::cin, line)) {
    lineCount++;

    /// Check if one minute has passed, if yes log the RAM usage
    auto currentTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::minutes>(currentTime - startTime).count();
    if (elapsedTime >= 1) {
      logRamUsage(_ramLogFile);
      startTime = currentTime;
    }

    /// progress indicator every 100k reads
    if (lineCount % 100000 == 0)
      std::cout << lineCount / 1000 << "k " << std::flush;

    if (line.empty() or startsWith(line, "@"))
      continue;

    /// split line into columns
    const auto fields = split(line);
    const auto& chromosome = fields[2];
    if (chromosome == "*" or fields[3] == "*" or fields[5] == "*")
      continue;

    validLines++;

    const auto cigar = parseCigar(fields[5]);
    const auto genomicRegion = regionRepresentation(chromosome , fields[3], cigar);
    const auto refseq = faidx(genomicRegion);

    const auto& readSeq = fields[9];
    const auto& readQual = fields[10];
    const std::string direction = (std::stoul(fields[1]) & 0x0010) ? "-" : "+";

    const auto candidateMutations = reconstructAlignment(cigar, readSeq, readQual, refseq, fields[3]);
    /// for all the candidate mutations found, update their stats if already found
    for (const auto& cluster : candidateMutations) {
      const auto refLength = cluster._refSeq.length();
      const auto varLength = cluster._readSeq.length();
      if (refLength >= 1 and varLength == 0)
        continue;

      if (refLength == 0 and varLength >= 1)
        continue;

      /// the unique key that represents a mutation, used for storing in main hash map
      const auto mutationKey = join("|", StringVector({chromosome, std::to_string(cluster._minPos), std::to_string(cluster._maxPos), cluster._refSeq, cluster._readSeq}));

      auto* foundStats = findMutationInCache(chromosome, cluster._minPos, mutationKey);
      if (foundStats == nullptr) {
        MutationStats stats;
        updateStats(stats, direction, cluster._meanQual);
        _mutations[convertChromosomeToNum(chromosome)][cluster._minPos][mutationKey] = stats;
        continue;
      }

      updateStats(*foundStats, direction, cluster._meanQual);
    }
  }

  std::cout << "Total lines: " << lineCount << std::endl
            << " - Valid lines: " << validLines << std::endl
            << " - Region cache hits: " << regionCacheHits << std::endl
            << " - Distinct mutations found: " << totalMutations() << std::endl;

  sam.close();
}

/// convert time into readable format
std::string hms(time_t seconds)
{
  auto sec = uint32_t (seconds % 60U);
  auto minutes = ((uint32_t) seconds - sec) / 60U;
  auto min = minutes % 60U;
  auto hours = (minutes - min) / 60U;

  std::stringstream ss;
  ss << hours
     << "h:" << std::setfill('0') << std::setw(2) << min
     << "m:" << std::setfill('0') << std::setw(2) << sec << "s";

  return ss.str();
}

/// main entry point of the program
int main(int argC, char* argV[])
{
  auto time0 = time(nullptr);

  IndelSeek indelSeek;
  indelSeek.parseArguments(argC, argV);
  indelSeek.printOptions();
  indelSeek.processFile();
  indelSeek.writeMutationsInVCF();

  auto time1 = time(nullptr);
  std::cout << "Done in " << hms(time1-time0) << std::endl;

  return 0;
}
