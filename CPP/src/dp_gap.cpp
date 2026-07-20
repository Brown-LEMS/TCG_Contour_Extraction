#include "dp_gap.hpp"

#include <cmath>
#include <limits>
#include <queue>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace tcg {
namespace {

struct Node {
  int x{};
  int y{};
  double theta{};
  double cost{};
  Node(int x_, int y_, double theta_, double c) : x(x_), y(y_), theta(theta_), cost(c) {}
};

struct NodeMinCmp {
  bool operator()(const Node* a, const Node* b) const { return a->cost > b->cost; }
};

double step_cost(const Node* p, const Node* q, const std::vector<double>& E_map,
                 const std::vector<double>& E_map_soft, const std::vector<double>& O_map, int h,
                 double sx, double sy, double contrast_th) {
  double wG = 0.65;
  double wS = 0.35;
  const int px = p->x;
  const int py = p->y;
  const int qx = q->x;
  const int qy = q->y;

  double fG = E_map_soft[static_cast<size_t>(qx * h + qy)];
  if (fG < contrast_th) return std::numeric_limits<double>::infinity();
  fG = fG * fG / (fG * fG + 0.01);

  double fZ = E_map[static_cast<size_t>(qx * h + qy)];
  if (fZ == 1.0) fZ = 0.0;

  const double dx = static_cast<double>(qx - px);
  const double dy = static_cast<double>(qy - py);
  const double dist = std::sqrt(dx * dx + dy * dy);

  double dq = std::cos(O_map[static_cast<size_t>(qx * h + qy)]) * dx +
              std::sin(O_map[static_cast<size_t>(qx * h + qy)]) * dy;
  dq = dq / dist;
  if (dq < 0) dq = -dq;
  const double dq_c = std::min(1.0, std::max(0.0, dq));
  const double fO = std::exp(-std::acos(dq_c) * std::acos(dq_c) / M_PI * 4.0 / M_PI * 4.0 / 2.0);

  const double cos_s =
      std::min(1.0, std::max(-1.0, (std::cos(p->theta) * dx + std::sin(p->theta) * dy) / dist));
  double fS = (1.0 / M_PI) * std::acos(cos_s);
  if (px == static_cast<int>(sx) && py == static_cast<int>(sy)) wS *= 2.0;

  return wG * (1.0 - fG * fO * (1.0 - fZ)) * dist + wS * fS * dist;
}

}  // namespace

DpGapResult dp_gap_cpt(const std::vector<double>& E_map, const std::vector<double>& E_map_soft,
                       const std::vector<double>& O_map, double start_x, double start_y, double start_theta,
                       int h, int w, double theta_th, double contrast_th) {
  DpGapResult out;
  const int N = h * w;
  out.cost.assign(static_cast<size_t>(N), std::numeric_limits<double>::infinity());
  out.back_p.assign(static_cast<size_t>(N), 0.0);
  out.len.assign(static_cast<size_t>(N), 0.0);

  std::vector<Node*> All_Nodes(static_cast<size_t>(N), nullptr);
  for (int x = 0; x < w; ++x) {
    for (int y = 0; y < h; ++y) {
      All_Nodes[static_cast<size_t>(x * h + y)] =
          new Node(x, y, O_map[static_cast<size_t>(x * h + y)], 0.0);
    }
  }

  const int sx = static_cast<int>(start_x);
  const int sy = static_cast<int>(start_y);
  const int sidx = sx * h + sy;
  out.cost[static_cast<size_t>(sidx)] = 0.0;
  out.len[static_cast<size_t>(sidx)] = 0.0;
  All_Nodes[static_cast<size_t>(sidx)]->theta = start_theta;

  using PriQ = std::priority_queue<Node*, std::vector<Node*>, NodeMinCmp>;
  PriQ L;
  std::vector<int> ind_L(static_cast<size_t>(N), 0);
  std::vector<int> ind_e(static_cast<size_t>(N), 0);

  L.push(All_Nodes[static_cast<size_t>(sidx)]);
  ind_L[static_cast<size_t>(sidx)] = 1;

  const double cos_th = std::cos(theta_th);

  while (!L.empty()) {
    Node* q_node = L.top();
    L.pop();
    const int qx = q_node->x;
    const int qy = q_node->y;
    const double qtheta = q_node->theta;
    ind_L[static_cast<size_t>(qx * h + qy)] = 0;
    ind_e[static_cast<size_t>(qx * h + qy)] = 1;

    std::vector<std::pair<int, int>> qualified_r;
    bool to_break = false;
    for (int dx = -1; dx <= 1; ++dx) {
      if (to_break) break;
      const int rx = qx + dx;
      if (rx < 0 || rx >= w) continue;
      for (int dy = -1; dy <= 1; ++dy) {
        if (to_break) break;
        const int ry = qy + dy;
        if (ry < 0 || ry >= h) continue;
        if (rx == qx && ry == qy) continue;
        if (ind_e[static_cast<size_t>(rx * h + ry)]) continue;

        const double dx0 = static_cast<double>(rx - sx);
        const double dy0 = static_cast<double>(ry - sy);
        const double n0 = std::sqrt(dx0 * dx0 + dy0 * dy0);
        if (n0 > 0) {
          const double cos_dir_diff0 = (std::cos(start_theta) * dx0 + std::sin(start_theta) * dy0) / n0;
          if (cos_dir_diff0 < cos_th) {
            ind_e[static_cast<size_t>(rx * h + ry)] = 1;
            continue;
          }
        }

        const double nd = std::sqrt(static_cast<double>(dx * dx + dy * dy));
        const double cos_dir_diff = (std::cos(qtheta) * dx + std::sin(qtheta) * dy) / nd;
        if (cos_dir_diff < 0.1) {
          ind_e[static_cast<size_t>(rx * h + ry)] = 1;
          continue;
        }

        if (E_map[static_cast<size_t>(rx * h + ry)] == 1.0 && !(rx == sx && ry == sy)) {
          qualified_r.clear();
          to_break = true;
        }
        qualified_r.emplace_back(dx, dy);
      }
    }

    for (const auto& dxy : qualified_r) {
      const int dx = dxy.first;
      const int dy = dxy.second;
      const int rx = qx + dx;
      const int ry = qy + dy;
      Node* r_node = All_Nodes[static_cast<size_t>(rx * h + ry)];
      const double cost_step =
          step_cost(q_node, r_node, E_map, E_map_soft, O_map, h, start_x, start_y, contrast_th);
      const double cost_temp = q_node->cost + cost_step;
      const double len_temp =
          out.len[static_cast<size_t>(qx * h + qy)] + std::sqrt(static_cast<double>(dx * dx + dy * dy));

      if (ind_L[static_cast<size_t>(rx * h + ry)] && cost_temp < r_node->cost) {
        out.cost[static_cast<size_t>(rx * h + ry)] = cost_temp;
        r_node->cost = cost_temp;
        r_node->theta = std::atan2(static_cast<double>(dy), static_cast<double>(dx));
        out.back_p[static_cast<size_t>(rx * h + ry)] = static_cast<double>(qx * h + qy);
        out.len[static_cast<size_t>(rx * h + ry)] = len_temp;
      } else if (!ind_L[static_cast<size_t>(rx * h + ry)]) {
        out.cost[static_cast<size_t>(rx * h + ry)] = cost_temp;
        r_node->cost = cost_temp;
        r_node->theta = std::atan2(static_cast<double>(dy), static_cast<double>(dx));
        out.back_p[static_cast<size_t>(rx * h + ry)] = static_cast<double>(qx * h + qy);
        out.len[static_cast<size_t>(rx * h + ry)] = len_temp;
        if (E_map[static_cast<size_t>(rx * h + ry)] == 1.0) continue;
        L.push(r_node);
        ind_L[static_cast<size_t>(rx * h + ry)] = 1;
      }
    }
  }

  for (Node* n : All_Nodes) delete n;
  out.cost[static_cast<size_t>(sidx)] = std::numeric_limits<double>::infinity();
  return out;
}

}  // namespace tcg
