#include "load_cem.hpp"

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

bool parse_paren_coords(const std::string& line, Edge& e) {
  auto l = line.find('(');
  auto r = line.find(')', l);
  if (l == std::string::npos || r == std::string::npos) return false;
  std::string inside = trim(line.substr(l + 1, r - l - 1));
  for (char& c : inside) {
    if (c == ',') c = ' ';
  }
  std::istringstream is(inside);
  if (!(is >> e.x >> e.y)) return false;
  std::string rest = line.substr(r + 1);
  std::istringstream rs(rest);
  if (!(rs >> e.dir >> e.conf >> e.d2f)) return false;
  return true;
}

}  // namespace

bool load_cem(const std::string& path, CemFile& out, std::string& error) {
  out = CemFile{};
  std::ifstream in(path);
  if (!in) {
    error = "Cannot open file: " + path;
    return false;
  }

  std::string line;
  if (!std::getline(in, line)) {
    error = "Empty file";
    return false;
  }
  line = trim(line);
  if (!starts_with(line, ".CEM v2.0")) {
    error = "Expected .CEM v2.0 header, got: " + line;
    return false;
  }

  enum class Section { None, Edgemap, Contours, Props };
  Section section = Section::None;
  int num_edges = 0;
  int num_contours = 0;
  size_t edges_read = 0;
  size_t props_read = 0;
  size_t next_contour_line = 0;

  while (std::getline(in, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#') continue;

    if (starts_with(line, "size=")) {
      int w = 0, h = 0;
      if (std::sscanf(line.c_str(), "size=[%d %d]", &w, &h) != 2) {
        error = "Bad size line: " + line;
        return false;
      }
      out.width = w;
      out.height = h;
      continue;
    }

    if (starts_with(line, "[Edgemap]")) {
      section = Section::Edgemap;
      continue;
    }
    if (starts_with(line, "[Contours]")) {
      section = Section::Contours;
      continue;
    }
    if (starts_with(line, "[Contour Properties]")) {
      section = Section::Props;
      continue;
    }

    if (section == Section::Edgemap) {
      if (starts_with(line, "count=")) {
        if (std::sscanf(line.c_str(), "count=%d", &num_edges) != 1) {
          error = "Bad edgemap count: " + line;
          return false;
        }
        out.edges.resize(static_cast<size_t>(num_edges));
        edges_read = 0;
        continue;
      }
      if (num_edges <= 0) continue;
      Edge e{};
      if (!parse_paren_coords(line, e)) {
        error = "Bad edgemap line: " + line;
        return false;
      }
      if (edges_read >= out.edges.size()) {
        error = "Too many edgemap rows";
        return false;
      }
      out.edges[edges_read++] = e;
      continue;
    }

    if (section == Section::Contours) {
      if (starts_with(line, "count=")) {
        if (std::sscanf(line.c_str(), "count=%d", &num_contours) != 1) {
          error = "Bad contour count: " + line;
          return false;
        }
        out.contours.assign(static_cast<size_t>(num_contours), {});
        out.contour_edge_idx.assign(static_cast<size_t>(num_contours), {});
        next_contour_line = 0;
        continue;
      }
      if (num_contours <= 0) continue;
      if (line.size() < 2 || line.front() != '[' || line.back() != ']') {
        error = "Bad contour index line: " + line;
        return false;
      }
      if (next_contour_line >= out.contours.size()) {
        error = "Too many contour lines";
        return false;
      }
      std::string inside = trim(line.substr(1, line.size() - 2));
      std::istringstream is(inside);
      ContourEdgeIndices ids;
      int id = 0;
      while (is >> id) ids.push_back(id);
      Contour chain;
      chain.reserve(ids.size());
      for (int eid : ids) {
        if (eid < 0 || static_cast<size_t>(eid) >= out.edges.size()) {
          error = "Contour references invalid edge id";
          return false;
        }
        chain.push_back(out.edges[static_cast<size_t>(eid)]);
      }
      out.contour_edge_idx[next_contour_line] = std::move(ids);
      out.contours[next_contour_line] = std::move(chain);
      ++next_contour_line;
      continue;
    }

    if (section == Section::Props) {
      if (num_contours <= 0) continue;
      if (line[0] == '#') continue;
      if (out.contour_props.empty()) {
        out.contour_props.resize(static_cast<size_t>(num_contours));
        props_read = 0;
      }
      if (props_read >= out.contour_props.size()) continue;
      std::istringstream ps(line);
      for (int j = 0; j < 8; ++j) {
        if (!(ps >> out.contour_props[props_read][j])) {
          error = "Bad contour property line: " + line;
          return false;
        }
      }
      ++props_read;
      continue;
    }
  }

  return true;
}

}  // namespace tcg
