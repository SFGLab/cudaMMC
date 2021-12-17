#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>
namespace fs = std::filesystem;

struct vector3 {
  float x, y, z;
};

constexpr std::string_view header_mmcif =
    "data_3dnome\n"
    "#\n"
    "_entry.id 3dnome\n"
    "#\n"
    "_audit_conform.dict_name       mmcif_pdbx.dic\n"
    "_audit_conform.dict_version    5.296\n"
    "_audit_conform.dict_location   "
    "http://mmcif.pdb.org/dictionaries/ascii/mmcif_pdbx.dic\n"
    "#\n"
    "loop_\n"
    "_atom_site.group_PDB\n"
    "_atom_site.id\n"
    "_atom_site.type_symbol\n"
    "_atom_site.label_atom_id\n"
    "_atom_site.label_alt_id\n"
    "_atom_site.label_comp_id\n"
    "_atom_site.label_asym_id\n"
    "_atom_site.label_entity_id\n"
    "_atom_site.label_seq_id\n"
    "_atom_site.pdbx_PDB_ins_code\n"
    "_atom_site.Cartn_x\n"
    "_atom_site.Cartn_y\n"
    "_atom_site.Cartn_z\n"
    "_atom_site.occupancy\n"
    "_atom_site.B_iso_or_equiv\n"
    "_atom_site.auth_asym_id";

std::string getWhiteSpaces(const int count) {
  std::stringstream output;
  for (int i = 0; i < count; ++i) {
    output << ' ';
  }
  return output.str();
}

std::string formatCoordinate(const float coordinate) {
  std::stringstream output;
  output.precision(3);

  if (coordinate > 0) {
    output << " " << coordinate << " ";
  } else {
    output << coordinate << " ";
  }

  const int length = output.str().length();

  if (length < 8) {
    output << getWhiteSpaces(9 - length);
  }

  return output.str();
}

std::vector<vector3> readData(const std::string &filename) {
  std::vector<vector3> result;
  std::ifstream file;
  file.exceptions(std::ifstream::badbit | std::ifstream::failbit);

  try {
    file.open(filename);
  } catch (const std::ifstream::failure &ex) {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(1);
  }

  std::string fullInput(std::istreambuf_iterator<char>{file}, {});
  std::stringstream inStream(fullInput);

  while (!inStream.eof()) {
    struct vector3 point;
    inStream >> point.x >> point.y >> point.z;
    result.push_back(point);
    inStream.ignore(1024, '\n');
  }

  return result;
}

void convert_to_pdb(std::string filename) {
  constexpr int extensionLength = 3;
  const std::string outputFileName = filename.replace(
      filename.length() - extensionLength, extensionLength, "pdb");
  std::ofstream outputFile(outputFileName);
  std::vector<vector3> coords = readData(filename);
  int i = 1;

  for (const auto &point : coords) {
    const auto x = formatCoordinate(point.x);
    const auto y = formatCoordinate(point.y);
    const auto z = formatCoordinate(-point.z);

    std::string strI = std::to_string(i);
    outputFile << "ATOM  " << getWhiteSpaces(5 - strI.length()) << strI
               << "   CA ALA A" << getWhiteSpaces(4 - strI.length()) << strI
               << "    ";
    outputFile << x << y << z << "  1.00 99.99\n";
    ++i;
  }

  outputFile.close();
}

// CONVERT TO MMCIF
void convert_to_mmcif(std::string filename) {
  auto atoms = std::vector<std::string>{std::string(header_mmcif)};
  std::string line;
  std::ifstream file;
  file.exceptions(std::ifstream::badbit | std::ifstream::failbit);

  // Open output stream
  constexpr int extensionLength = 3;
  const std::string outputFileName = filename.replace(
      filename.length() - extensionLength, extensionLength, "pdb");
  std::ofstream output_file(outputFileName);

  try {
    file.open(filename);
  } catch (const std::ifstream::failure &ex) {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(1);
  }
  while (getline(file, line)) {
    std::istringstream iss(line);
    std::vector<std::string> values{std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>()};
    auto i = 1;
    std::string num;
    if (values[4] != "A") {
      i = 0;
      const auto a_loc = values[4].find("A");
      num = values[4].substr(a_loc);
    } else {
      num = values[5];
    }
    output_file << "ATOM " + values[1] + " C " + values[2] + " . " + values[3] +
                       " A  1 " + num + " ? " + values[5 + i] + " " +
                       values[6 + i] + " " + values[7 + i] + " " +
                       values[8 + i] + " " + values[9 + i] + " A\n";
  }
  // OFSTREAM IS RAII SO NO NEED TO CLOSE
}

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cout << "buongiorno ragazzi\n"
              << " Please provide 2 arguments\n"
              << "     1) in_directory string\n"
              << "     2) smooth num" << std::endl;
    return -1;
  }
  const auto in_directory = std::string(argv[1]);
  const auto smooth_num = std::string(argv[2]);

  std::vector<std::string> hcms;
  std::vector<std::string> txts;
  for (const auto &entry : fs::directory_iterator(in_directory)) {
    const auto tmp_path = std::string(entry.path());
    if (tmp_path.substr(tmp_path.size() - 3).compare(std::string("hcm")) == 0) {
      hcms.push_back(tmp_path);
    }
  }

  // Generate smooth files
  for (const auto &hcm : hcms) {
    const std::string command =
        "./cuMMC -a smooth -i " + hcm + " -r " + smooth_num;
    std::system(command.c_str());
  }

  // Get all smooth files
  for (const auto &entry : fs::directory_iterator(in_directory)) {
    const auto tmp_path = std::string(entry.path());
    if (tmp_path.substr(tmp_path.size() - 10)
            .compare(std::string("smooth.txt")) == 0) {
      txts.push_back(tmp_path);
    }
  }

  // convert
  for (const auto &txt : txts) {
    convert_to_pdb(txt);
    convert_to_mmcif(txt.substr(0, txt.size() - 3) + "pdb");
  }

  return 0;
}
