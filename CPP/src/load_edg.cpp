#include "load_edg.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

namespace tcg {

namespace {

std::string trim(std::string s) {
  auto a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return {};
  auto b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

bool starts_with(const std::string& s, const char* p) {
  size_t n = std::char_traits<char>::length(p);
  return s.size() >= n && s.compare(0, n, p) == 0;
}

/// Parse "[a, b]" or "(a, b)" starting at position `pos`; advances `pos` past closing bracket.
bool parse_pair_brackets(const std::string& line, size_t& pos, double& a, double& b) {
  size_t open = line.find_first_of("[(", pos);
  if (open == std::string::npos) return false;
  char close = (line[open] == '[') ? ']' : ')';
  size_t shut = line.find(close, open + 1);
  if (shut == std::string::npos) return false;
  std::string inside = trim(line.substr(open + 1, shut - open - 1));
  for (char& c : inside)
    if (c == ',') c = ' ';
  std::istringstream is(inside);
  if (!(is >> a >> b)) return false;
  pos = shut + 1;
  return true;
}

bool parse_pair_brackets_int(const std::string& line, size_t& pos, int& a, int& b) {
  double da, db;
  if (!parse_pair_brackets(line, pos, da, db)) return false;
  a = static_cast<int>(da);
  b = static_cast<int>(db);
  return true;
}

}  // namespace

bool load_edg(const std::string& path, EdgFile& out, std::string& error) {
  out = EdgFile{};
  std::ifstream in(path);
  if (!in) {
    error = "Cannot open file: " + path;
    return false;
  }

  int ver = 1;
  int width = 0, height = 0;
  int num_edges = 0;
  bool have_count = false;
  std::vector<Edge> edges;
  size_t next_edge = 0;

  std::string line;
  while (std::getline(in, line)) {
    line = trim(line);
    if (line.empty()) continue;

    if (starts_with(line, "# EDGE_MAP v2.0") || starts_with(line, "# EDGE_MAP v3.0")) {
      ver = 2;
      continue;
    }
    if (line[0] == '#') continue;

    if (starts_with(line, "WIDTH=")) {
      if (std::sscanf(line.c_str(), "WIDTH=%d", &width) != 1 &&
          std::sscanf(line.c_str(), " WIDTH=%d", &width) != 1) {
        error = "Bad WIDTH line: " + line;
        return false;
      }
      continue;
    }
    if (starts_with(line, "HEIGHT=")) {
      if (std::sscanf(line.c_str(), "HEIGHT=%d", &height) != 1 &&
          std::sscanf(line.c_str(), " HEIGHT=%d", &height) != 1) {
        error = "Bad HEIGHT line: " + line;
        return false;
      }
      continue;
    }
    if (starts_with(line, "EDGE_COUNT=")) {
      if (std::sscanf(line.c_str(), "EDGE_COUNT=%d", &num_edges) != 1 &&
          std::sscanf(line.c_str(), " EDGE_COUNT=%d", &num_edges) != 1) {
        error = "Bad EDGE_COUNT line: " + line;
        return false;
      }
      edges.resize(num_edges);
      next_edge = 0;
      have_count = true;
      out.edgemap.assign(static_cast<size_t>(height) * static_cast<size_t>(width), 0.0);
      out.thetamap.assign(static_cast<size_t>(height) * static_cast<size_t>(width), 0.0);
      continue;
    }

    if (!have_count) continue;

    size_t bracket = line.find('[');
    if (bracket == std::string::npos) continue;

    size_t p = bracket;
    int ix = 0, iy = 0;
    if (!parse_pair_brackets_int(line, p, ix, iy)) continue;

    double idir = 0, iconf = 0;
    {
      std::istringstream is(line.substr(p));
      if (!(is >> idir >> iconf)) continue;
      p += static_cast<size_t>(is.tellg());
    }

    double x = 0, y = 0;
    if (!parse_pair_brackets(line, p, x, y)) continue;

    Edge e{};
    e.x = x;
    e.y = y;
    if (ver == 1) {
      if (std::sscanf(line.c_str() + p, "%lf %lf", &e.dir, &e.conf) < 2) continue;
      if (e.conf == 0) e.conf = 2;
      e.d2f = 0;
    } else {
      int n = std::sscanf(line.c_str() + p, "%lf %lf %lf", &e.dir, &e.conf, &e.d2f);
      if (n < 2) continue;
      if (n < 3) e.d2f = 0;
      if (e.conf == 0) e.conf = 2;
    }

    if (next_edge >= edges.size()) {
      error = "Too many EDGE lines";
      return false;
    }
    edges[next_edge] = e;

    if (width > 0 && height > 0 && !out.edgemap.empty()) {
      int row = static_cast<int>(std::lround(e.y + 1.0)) - 1;
      int col = static_cast<int>(std::lround(e.x + 1.0)) - 1;
      if (row >= 0 && row < height && col >= 0 && col < width) {
        size_t idx = static_cast<size_t>(row) * static_cast<size_t>(width) + static_cast<size_t>(col);
        out.edgemap[idx] = e.conf;
        out.thetamap[idx] = e.dir;
      }
    }
    ++next_edge;
  }

  if (have_count && next_edge != edges.size()) {
    error = "EDGE_COUNT=" + std::to_string(static_cast<int>(edges.size())) + " but read " +
            std::to_string(next_edge) + " edge lines";
    return false;
  }

  out.version = ver;
  out.width = width;
  out.height = height;
  out.edges = std::move(edges);
  return true;
}

}  // namespace tcg
