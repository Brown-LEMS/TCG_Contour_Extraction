#include "matlab_bwmorph.hpp"
#include "matlab_lutskel.hpp"

namespace tcg {
namespace {

// MATLAB bwlookup 3x3 bit weights (column-major):
//   1   8   64
//   2  16  128
//   4  32  256
inline int bwlookup_index3_fg(const uchar* row0, const uchar* row1, const uchar* row2, int x, int w) {
  // Center is known foreground (bit 16).
  auto at = [&](const uchar* row, int xx) -> int {
    if (!row || xx < 0 || xx >= w) return 0;
    return row[xx] ? 1 : 0;
  };
  return at(row0, x - 1) * 1 + at(row1, x - 1) * 2 + at(row2, x - 1) * 4 + at(row0, x) * 8 + 16 +
         at(row2, x) * 32 + at(row0, x + 1) * 64 + at(row1, x + 1) * 128 + at(row2, x + 1) * 256;
}

// Apply one lutskel LUT. Background stays 0 (LUT encodes center & ~hitmiss).
// Returns whether any foreground pixel was deleted.
bool apply_lutskel(const cv::Mat& in, cv::Mat& out, const uint8_t lut[512]) {
  const int h = in.rows;
  const int w = in.cols;
  out.setTo(0);
  bool changed = false;
  for (int y = 0; y < h; ++y) {
    const uchar* r0 = (y > 0) ? in.ptr<uchar>(y - 1) : nullptr;
    const uchar* r1 = in.ptr<uchar>(y);
    const uchar* r2 = (y + 1 < h) ? in.ptr<uchar>(y + 1) : nullptr;
    uchar* o = out.ptr<uchar>(y);
    for (int x = 0; x < w; ++x) {
      if (!r1[x]) continue;  // background stays 0
      const uchar v = lut[bwlookup_index3_fg(r0, r1, r2, x, w)] ? 1 : 0;
      o[x] = v;
      if (!v) changed = true;
    }
  }
  return changed;
}

}  // namespace

cv::Mat morph_skel_matlab(const cv::Mat& bin) {
  // Port of MATLAB: bwmorph(BW,'skel',Inf) / images.internal.algbwmorph 'skeleton'.
  const int h = bin.rows;
  const int w = bin.cols;
  cv::Mat a(h, w, CV_8UC1);
  cv::Mat b(h, w, CV_8UC1);
  for (int y = 0; y < h; ++y) {
    const uchar* s = bin.ptr<uchar>(y);
    uchar* d = a.ptr<uchar>(y);
    for (int x = 0; x < w; ++x) d[x] = s[x] ? 1 : 0;
  }

  cv::Mat* cur = &a;
  cv::Mat* nxt = &b;
  for (;;) {
    bool any = false;
    for (int i = 0; i < matlab_bwmorph::kLutSkelCount; ++i) {
      any = apply_lutskel(*cur, *nxt, matlab_bwmorph::kLutSkel[i]) || any;
      cv::Mat* tmp = cur;
      cur = nxt;
      nxt = tmp;
    }
    if (!any) break;
  }
  return cur->clone();
}

}  // namespace tcg
