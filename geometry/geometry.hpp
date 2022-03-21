#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
namespace Geometry {

class Vector {
 private:
  std::pair<int64_t, int64_t> vector_coord_;

 public:
  Vector();
  Vector(std::pair<int64_t, int64_t> begin, std::pair<int64_t, int64_t> end);
  Vector(std::pair<int64_t, int64_t> coord);
  double GetLength() const;
  std::pair<int64_t, int64_t> GetCoord() const;
  Vector& operator+=(Vector& added);
  Vector operator+(Vector& add);
  int64_t operator*(const Vector& prod) const;

  int64_t operator^(const Vector& sec) const;
  Vector operator-();
  Vector& operator-=(Vector& subtr);
  Vector operator-(Vector& subtr);
  ~Vector() = default;
};
class Point;
class Segment;
class IShape {
 public:
  virtual IShape& Move(const Vector& vt) = 0;
  virtual bool ContainsPoint(const Point& pt) const = 0;
  virtual bool CrossesSegment(const Segment&) const = 0;
  virtual IShape* Clone() const = 0;
  virtual std::string ToString() = 0;
  virtual ~IShape() = default;
};
class Point : public IShape {
 public:
  int64_t x_coord = 0;
  int64_t y_coord = 0;
  Point();
  Point(int64_t x, int64_t y);
  virtual IShape& Move(const Vector& vt) override;
  virtual bool ContainsPoint(const Point& pt) const override;
  virtual bool CrossesSegment(const Segment& sg) const override;
  Point* Clone() const override;
  virtual std::string ToString() override;
  virtual ~Point() override = default;
};
class Segment : public IShape {
 private:
  std::pair<int64_t, int64_t> beg_;
  std::pair<int64_t, int64_t> end_;

 public:
  Segment();
  Segment(const std::pair<int64_t, int64_t>& beg,
          const std::pair<int64_t, int64_t>& end);
  Segment(const Point& pt1, const Point& pt2);
  double DistToSeg(std::pair<int64_t, int64_t>& pt) const;
  virtual IShape& Move(const Vector& vt) override;
  virtual bool ContainsPoint(const Point& pt) const override;
  virtual bool CrossesSegment(const Segment& sg) const override;
  Segment* Clone() const override;
  virtual std::string ToString() override;
  std::pair<int64_t, int64_t> GetBeg() const;
  std::pair<int64_t, int64_t> GetEnd() const;
  virtual ~Segment() override = default;
};
class Ray : public IShape {
 private:
  std::pair<int64_t, int64_t> beg_;
  std::pair<int64_t, int64_t> pt_on_ray_;

 public:
  Ray();
  Ray(const std::pair<int64_t, int64_t>& beg,
      const std::pair<int64_t, int64_t>& pt_on_ray);
  Ray(const Point& beg, const Point& end);
  double DistToRay(std::pair<int64_t, int64_t>& pt) const;
  virtual IShape& Move(const Vector& vt) override;
  virtual bool ContainsPoint(const Point& pt) const override;
  virtual bool CrossesSegment(const Segment& sg) const override;
  Ray* Clone() const override;
  virtual std::string ToString() override;
  virtual ~Ray() override = default;
};
class Line : public IShape {
 private:
  int64_t A_ = 0;
  int64_t B_ = 0;
  int64_t C_ = 0;
  Point pt1_;
  Point pt2_;
  // TODO: two points
 public:
  Line(int64_t a, int64_t b, int64_t c);
  Line(std::pair<int64_t, int64_t> pt1, std::pair<int64_t, int64_t> pt2);
  Line(const Point& pt1, const Point& pt2);
  Line();
  Vector DirectionVector() const;
  void Intersection(Line& sec) const;
  double DistanceBetweenParall(Line& sec) const;
  bool IsParall(const Line& sec) const;
  double DistToLine(Point& pt) const;
  virtual IShape& Move(const Vector& vt) override;
  virtual bool ContainsPoint(const Point& pt) const override;
  virtual bool CrossesSegment(const Segment& sg) const override;
  Line* Clone() const override;
  virtual std::string ToString() override;
  virtual ~Line() override = default;
};
class Circle : public IShape {
 private:
  Point center_;
  int64_t radius_;

 public:
  Circle();
  Circle(Point center, int64_t rad);
  virtual IShape& Move(const Vector& vt) override;
  virtual bool ContainsPoint(const Point& pt) const override;
  virtual bool CrossesSegment(const Segment& sg) const override;
  Circle* Clone() const override;
  virtual std::string ToString() override;
  virtual ~Circle() override = default;
};
class Polygon : public IShape {
 private:
  uint64_t vertexes_ = 0;
  std::vector<Point> coord_;

 public:
  Polygon();
  Polygon(uint64_t vert, std::vector<Point> coord);
  Polygon(std::vector<Point> coord);
  virtual IShape& Move(const Vector& vt) override;
  virtual bool ContainsPoint(const Point& pt) const override;
  virtual bool CrossesSegment(const Segment& sg) const override;
  Polygon* Clone() const override;
  virtual std::string ToString() override;
  virtual ~Polygon() override = default;
};
Vector::Vector() = default;
Vector::Vector(std::pair<int64_t, int64_t> begin,
               std::pair<int64_t, int64_t> end) {
  vector_coord_.first = end.first - begin.first;
  vector_coord_.second = end.second - begin.second;
}
Vector::Vector(std::pair<int64_t, int64_t> coord) { vector_coord_ = coord; }
double Vector::GetLength() const {
  double ans = sqrt(vector_coord_.first * vector_coord_.first +
                    vector_coord_.second * vector_coord_.second);
  return ans;
}
std::pair<int64_t, int64_t> Vector::GetCoord() const { return vector_coord_; }
Vector& Vector::operator+=(Vector& added) {
  vector_coord_.first += added.vector_coord_.first;
  vector_coord_.second += added.vector_coord_.second;
  return *this;
}
Vector Vector::operator+(Vector& add) {
  Vector copy = *this;
  copy += add;
  return copy;
}
int64_t Vector::operator*(const Vector& prod) const {
  return (vector_coord_.first * prod.vector_coord_.first +
          vector_coord_.second * prod.vector_coord_.second);
}
int64_t Vector::operator^(const Vector& sec) const {
  return (vector_coord_.first * sec.vector_coord_.second -
          vector_coord_.second * sec.vector_coord_.first);
}
Vector Vector::operator-() {
  std::pair<int64_t, int64_t> beg = {0, 0};
  std::pair<int64_t, int64_t> end = {-vector_coord_.first,
                                     -vector_coord_.second};
  return Vector(beg, end);
}
Vector& Vector::operator-=(Vector& subtr) {
  vector_coord_.first -= subtr.vector_coord_.first;
  vector_coord_.second -= subtr.vector_coord_.second;
  return *this;
}
Vector Vector::operator-(Vector& subtr) {
  Vector copy = *this;
  copy -= subtr;
  return copy;
}
Point::Point() = default;
Point::Point(int64_t x, int64_t y) {
  x_coord = x;
  y_coord = y;
}
IShape& Point::Move(const Vector& vt) {
  x_coord += vt.GetCoord().first;
  y_coord += vt.GetCoord().second;
  return *this;
}
bool Point::ContainsPoint(const Point& pt) const {
  return (pt.x_coord == x_coord && pt.y_coord == y_coord);
}
bool Point::CrossesSegment(const Segment& sg) const {
  return sg.ContainsPoint(*this);
}
Point* Point::Clone() const {
  Point* copy = new Point(x_coord, y_coord);
  return copy;
}
std::string Point::ToString() {
  std::string s1 =
      "Point(" + std::to_string(x_coord) + ", " + std::to_string(y_coord) + ")";
  return s1;
}
Vector operator-(const Point& red, const Point& substr) {
  std::pair<int64_t, int64_t> a = {red.x_coord, red.y_coord};
  std::pair<int64_t, int64_t> b = {substr.x_coord, substr.y_coord};
  Vector v1(b, a);
  return v1;
}
Line::Line() = default;
Line::Line(int64_t a, int64_t b, int64_t c) {
  this->A_ = a;
  this->B_ = b;
  this->C_ = c;
}
Line::Line(std::pair<int64_t, int64_t> pt1, std::pair<int64_t, int64_t> pt2) {
  this->A_ = pt2.second - pt1.second;
  this->B_ = pt1.first - pt2.first;
  this->C_ = pt1.second * (pt2.first - pt1.first) -
             pt1.first * (pt2.second - pt1.second);
  Point a(pt1.first, pt1.second);
  Point b(pt2.first, pt1.second);
}
Line::Line(const Point& pt1, const Point& pt2) {
  this->A_ = pt2.y_coord - pt1.y_coord;
  this->B_ = pt1.x_coord - pt2.x_coord;
  this->C_ = pt1.y_coord * (pt2.x_coord - pt1.x_coord) -
             pt1.x_coord * (pt2.y_coord - pt1.y_coord);
  this->pt1_ = pt1;
  this->pt2_ = pt2;
}
Vector Line::DirectionVector() const {
  std::pair<int64_t, int64_t> beg = {0, 0};
  std::pair<int64_t, int64_t> end = {B_, -A_};
  return Vector(beg, end);
}
void Line::Intersection(Line& sec) const {
  std::cout << (((-1) * C_ * sec.B_ + sec.C_ * B_) /
                (A_ * sec.B_ - sec.A_ * B_))
            << (((-1) * A_ * sec.C_ + sec.A_ * C_) /
                (A_ * sec.B_ - sec.A_ * B_))
            << std::endl;
}
double Line::DistanceBetweenParall(Line& sec) const {
  return std::fabs((C_ - (A_ / sec.A_) * sec.C_) / sqrt(A_ * A_ + B_ * B_));
}
bool Line::IsParall(const Line& sec) const {
  return ((A_ == sec.A_ && B_ == sec.B_) || (A_ * sec.B_ == sec.A_ * B_));
}
double Line::DistToLine(Point& pt) const {
  int64_t ans = std::abs(A_ * pt.x_coord + B_ * pt.y_coord + C_);
  return (static_cast<double>(ans) / sqrt(A_ * A_ + B_ * B_));
}
IShape& Line::Move(const Vector& vt) {
  C_ = C_ - A_ * vt.GetCoord().first - B_ * vt.GetCoord().second;
  return *this;
}
bool Line::ContainsPoint(const Point& pt) const {
  int64_t ans = A_ * pt.x_coord + B_ * pt.y_coord + C_;
  return ans == 0;
}
bool Line::CrossesSegment(const Segment& sg) const {
  Vector sg_v(sg.GetBeg(), sg.GetEnd());
  Vector dir = DirectionVector();
  std::pair<int64_t, int64_t> a(pt1_.x_coord, pt1_.y_coord);
  Vector v1(a, sg.GetBeg());
  Vector v2(a, sg.GetEnd());
  return (((dir ^ v1) * (dir ^ v2)) <= 0);
}
Line* Line::Clone() const {
  Line* copy = new Line(A_, B_, C_);
  return copy;
}
std::string Line::ToString() {
  std::string s1 = "Line(" + std::to_string(A_) + ", " + std::to_string(B_) +
                   ", " + std::to_string(C_) + ")";
  return s1;
}
Ray::Ray() = default;
Ray::Ray(const std::pair<int64_t, int64_t>& beg,
         const std::pair<int64_t, int64_t>& pt_on_ray) {
  this->beg_ = beg;
  this->pt_on_ray_ = pt_on_ray;
}
Ray::Ray(const Point& beg, const Point& end) {
  beg_.first = beg.x_coord;
  beg_.second = beg.y_coord;
  pt_on_ray_.first = end.x_coord;
  pt_on_ray_.second = end.y_coord;
}
bool Ray::ContainsPoint(const Point& pt) const {
  std::pair<int64_t, int64_t> a(pt.x_coord, pt.y_coord);
  Vector v1(beg_, a);
  Vector v2(beg_, pt_on_ray_);
  Line l1(beg_, pt_on_ray_);
  return (v1 * v2 >= 0 && l1.ContainsPoint(pt));
}
double Ray::DistToRay(std::pair<int64_t, int64_t>& pt) const {
  Vector v1(beg_, pt_on_ray_);
  Vector v2(beg_, pt);
  Line l1 = Line(beg_, pt_on_ray_);
  int64_t sc_mul = v1 * v2;
  Point a(pt.first, pt.second);
  return sc_mul >= 0 ? l1.DistToLine(a) : v2.GetLength();
}
IShape& Ray::Move(const Vector& vt) {
  Point pt(beg_.first, beg_.second);
  pt.Move(vt);
  Point end(pt_on_ray_.first, pt_on_ray_.second);
  end.Move(vt);
  beg_.first = pt.x_coord;
  beg_.second = pt.y_coord;
  pt_on_ray_.first = end.x_coord;
  pt_on_ray_.second = end.y_coord;
  return *this;
}
template <typename T>
int Sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
bool Ray::CrossesSegment(const Segment& sg) const {
  Vector v1 = Vector(beg_, pt_on_ray_);
  Vector v2 = Vector(sg.GetBeg(), sg.GetEnd());
  Vector v3 = Vector(sg.GetBeg(), beg_);
  if (Sgn(v1 ^ v2) == 0) {
    Point sg_beg(sg.GetBeg().first, sg.GetBeg().second);
    Point sg_pt(sg.GetEnd().first, sg.GetEnd().second);
    return ContainsPoint(sg_beg) || ContainsPoint(sg_pt);
  }
  int64_t rs = v1 ^ v2;
  int64_t rsr = v1 ^ v3;
  int64_t srs = v2 ^ v3;
  Sgn(rs) < 0 ? (rs *= -1, rsr *= -1, srs *= -1) : 0;
  return rsr <= rs && rsr >= 0 && srs >= 0;
}
Ray* Ray::Clone() const {
  Ray* clone = new Ray(beg_, pt_on_ray_);
  return clone;
}
std::string Ray::ToString() {
  Point pt(pt_on_ray_.first, pt_on_ray_.second);
  Point beg(beg_.first, beg_.second);
  Vector tmp = pt - beg;
  std::string s1 = "Ray(Point(" + std::to_string(beg_.first) + ", " +
                   std::to_string(beg_.second) + "), " + "Vector(" +
                   std::to_string(tmp.GetCoord().first) + ", " +
                   std::to_string(tmp.GetCoord().second) + "))";
  return s1;
}
Segment::Segment() = default;
Segment::Segment(const std::pair<int64_t, int64_t>& beg,
                 const std::pair<int64_t, int64_t>& end)
    : beg_(beg), end_(end) {}
Segment::Segment(const Point& pt1, const Point& pt2) {
  std::pair<int64_t, int64_t> a(pt1.x_coord, pt1.y_coord);
  std::pair<int64_t, int64_t> b(pt2.x_coord, pt2.y_coord);
  beg_ = a;
  end_ = b;
}
bool Segment::ContainsPoint(const Point& pt) const {
  Vector vect_ac = Vector({pt.x_coord - beg_.first, pt.y_coord - beg_.second});
  Vector vect_ab = Vector(beg_, end_);
  return (vect_ab ^ vect_ac) == 0 &&
         pt.x_coord <= std::max(beg_.first, end_.first) &&
         pt.x_coord >= std::min(beg_.first, end_.first) &&
         pt.y_coord <= std::max(beg_.second, end_.second) &&
         pt.y_coord >= std::min(beg_.second, end_.second);
}
std::pair<int64_t, int64_t> Segment::GetBeg() const { return beg_; }
std::pair<int64_t, int64_t> Segment::GetEnd() const { return end_; }
double Segment::DistToSeg(std::pair<int64_t, int64_t>& pt) const {
  Vector v1(beg_, pt);
  Vector v2(beg_, end_);
  Vector v3(end_, pt);
  int64_t v1_v2 = v1 * v2;
  int64_t v3_v4 = v3 * (-v2);
  Line l1(beg_, end_);
  Point a(pt.first, pt.second);
  return (v1_v2 >= 0 && v3_v4 >= 0) ? l1.DistToLine(a)
                                    : std::min(v1.GetLength(), v3.GetLength());
}
IShape& Segment::Move(const Vector& vt) {
  beg_.first += vt.GetCoord().first;
  beg_.second += vt.GetCoord().second;
  end_.first += vt.GetCoord().first;
  end_.second += vt.GetCoord().second;
  return *this;
}
bool Segment::CrossesSegment(const Segment& sg) const {
  Vector ab_vect(beg_, end_);
  Vector ad_vect(beg_, sg.end_);
  Vector ac_vect(beg_, sg.beg_);
  Vector cd_vect(sg.beg_, sg.end_);
  Vector ca_vect(sg.beg_, beg_);
  Vector cb_vect(sg.beg_, end_);
  Vector da_vect(sg.end_, beg_);
  Vector db_vect(sg.end_, end_);
  Vector bd_vect(end_, sg.end_);
  Vector bc_vect(end_, sg.beg_);
  return (((ab_vect ^ ac_vect) * (ab_vect ^ ad_vect)) < 0 &&
          ((cd_vect ^ ca_vect) * (cd_vect ^ cb_vect)) < 0) ||
         ((ab_vect ^ ac_vect) == 0 && (ca_vect) * (cb_vect) <= 0) ||
         ((cd_vect ^ ca_vect) == 0 && (ac_vect) * (ad_vect) <= 0) ||
         ((cd_vect ^ cb_vect) == 0 && (bc_vect) * (bd_vect) <= 0) ||
         ((ab_vect ^ ad_vect) == 0 && (da_vect) * (db_vect) <= 0);
}
Segment* Segment::Clone() const {
  Segment* sg = new Segment(beg_, end_);
  return sg;
}
std::string Segment::ToString() {
  std::string s1 = "Segment(Point(" + std::to_string(beg_.first) + ", " +
                   std::to_string(beg_.second) + "), " + "Point(" +
                   std::to_string(end_.first) + ", " +
                   std::to_string(end_.second) + "))";
  return s1;
}
Polygon::Polygon() = default;
Polygon::Polygon(uint64_t vert, std::vector<Point> coord) {
  vertexes_ = vert;
  for (size_t i = 0; i < coord.size(); ++i) {
    coord_.push_back(coord[i]);
  }
}
Polygon::Polygon(std::vector<Point> coord) {
  for (size_t i = 0; i < coord.size(); ++i) {
    coord_.push_back(coord[i]);
  }
}
bool Polygon::ContainsPoint(const Point& pt) const {
  for (size_t i = 0; i < coord_.size(); ++i) {
    uint64_t last = (i + 1) % coord_.size();
    Segment sg(coord_[i], coord_[last]);
    if (sg.ContainsPoint(pt)) {
      return true;
    }
  }
  Point far(1e5 + 1, pt.y_coord + 1);
  size_t counter = 0;
  Segment s1(pt, far);
  for (size_t i = 0; i < coord_.size(); ++i) {
    uint64_t last = (i + 1) % coord_.size();
    Segment s2(coord_[i], coord_[last]);
    if (s1.CrossesSegment(s2)) {
      ++counter;
    }
  }
  return ((counter % 2) == 1);
}
bool Polygon::CrossesSegment(const Segment& sg) const {
  for (size_t i = 0; i < coord_.size(); ++i) {
    uint64_t last = (i + 1) % coord_.size();
    ;
    Segment s2(coord_[i], coord_[last]);
    Segment s3(coord_[coord_.size() - 1], coord_[0]);
    if (s2.CrossesSegment(sg)) {
      return true;
    }
    if (s3.CrossesSegment(sg)) {
      return true;
    }
  }
  return false;
}

Polygon* Polygon::Clone() const {
  Polygon* copy = new Polygon(vertexes_, coord_);
  return copy;
}
std::string Polygon::ToString() {
  std::string s1 = "Polygon(";
  for (size_t i = 0; i < coord_.size(); ++i) {
    s1 += "Point(" + std::to_string(coord_[i].x_coord) + ", " +
          std::to_string(coord_[i].y_coord) + ")";
    if (i != coord_.size() - 1) {
      s1 += ", ";
    }
  }
  s1 += ")";
  return s1;
}
IShape& Polygon::Move(const Vector& vt) {
  for (size_t i = 0; i < coord_.size(); ++i) {
    coord_[i].x_coord += vt.GetCoord().first;
    coord_[i].y_coord += vt.GetCoord().second;
  }
  return *this;
}
Circle::Circle() = default;
Circle::Circle(Point center, int64_t rad) {
  center_ = center;
  radius_ = rad;
}
IShape& Circle::Move(const Vector& vt) {
  center_.x_coord += vt.GetCoord().first;
  center_.y_coord += vt.GetCoord().second;
  return *this;
}
bool Circle::ContainsPoint(const Point& pt) const {
  if (std::abs(pt.x_coord - center_.x_coord) > radius_) {
    return false;
  }
  if (std::abs(pt.y_coord - center_.y_coord) > radius_) {
    return false;
  }
  if (std::abs(pt.x_coord - center_.x_coord) +
          std::abs(pt.y_coord - center_.y_coord) <=
      radius_) {
    return true;
  }
  int32_t tmp_x = std::abs(pt.x_coord - center_.x_coord);
  int32_t tmp_y = std::abs(pt.y_coord - center_.y_coord);
  return (tmp_x * tmp_x + tmp_y * tmp_y <= radius_ * radius_);
}
bool Circle::CrossesSegment(const Segment& sg) const {
  Point a(sg.GetBeg().first, sg.GetBeg().second);
  Point b(sg.GetEnd().first, sg.GetEnd().second);
  Vector v1(sg.GetBeg(), sg.GetEnd());
  std::pair<int64_t, int64_t> s1(center_.x_coord, center_.y_coord);
  int64_t dist = sg.DistToSeg(s1);
  Vector v_c1(s1, sg.GetBeg());
  Vector v_c2(s1, sg.GetEnd());
  return dist <= radius_ &&
         (v_c1.GetLength() >= radius_ || v_c2.GetLength() >= radius_);
}
std::string Circle::ToString() {
  std::string s1 = "Circle(Point(" + std::to_string(center_.x_coord) + ", " +
                   std::to_string(center_.y_coord) + "), " +
                   std::to_string(radius_) + ")";
  return s1;
}
double TriangleSquare(Vector& fir, Vector& sec) {
  return static_cast<double>(std::abs(fir ^ sec)) / 2;
}
Circle* Circle::Clone() const {
  Circle* circ = new Circle(center_, radius_);
  return circ;
}
}  // namespace Geometry
template <class SmartPtrT>
void Delete(const SmartPtrT& dlt) {}

template <class T>
void Delete(T* ptr) {
  delete ptr;
}

void CheckFunctions(const Geometry::IShape* shape_ptr,
                    const Geometry::Point& point_a,
                    const Geometry::Point& point_b) {
  std::cout << "Given shape "
            << (shape_ptr->ContainsPoint(point_a) ? "contains"
                                                  : "does not contain")
            << " point A\n";

  const auto kSegmentAb = Geometry::Segment(point_a, point_b);
  std::cout << "Given shape "
            << (shape_ptr->CrossesSegment(kSegmentAb) ? "crosses"
                                                      : "does not cross")
            << " segment AB\n";

  const auto kVectorAb = point_b - point_a;
  const auto kClonedShapePtr =
      shape_ptr->Clone();  // may return either raw or smart pointer
  std::cout << kClonedShapePtr->Move(kVectorAb).ToString();

  Delete(kClonedShapePtr);  // raw pointer compatibility
}