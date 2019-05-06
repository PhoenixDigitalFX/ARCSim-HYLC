#include "HermiteSpline2D.hpp"

HermiteSpline2D::Mat4x4 HermiteSpline2D::M = [](){
    // M = np.array([
    //         [2,-2,1,1],
    //         [-3,3,-2,-1],
    //         [0,0,1,0],
    //         [1,0,0,0]
    //     ])
    return Mat4x4 {Vec4{2,-3,0,1},Vec4{-2,3,0,0},Vec4{1,-2,1,0},Vec4{1,-1,0,0}};
}();
HermiteSpline2D::Mat4x4 HermiteSpline2D::MT = [](){
    return (Mat4x4 {Vec4{2,-3,0,1},Vec4{-2,3,0,0},Vec4{1,-2,1,0},Vec4{1,-1,0,0}}).t();
}();


HermiteSpline2D::Mat4x4 HermiteSpline2D::MextL = [](){
    // [0,0,0,0],
    // [0,0,0,0],
    // [0,0,lin,0],
    // [1,0,0,0]
    return (Mat4x4 {Vec4{0,0,0,1},Vec4{0,0,0,0},Vec4{0,0,1,0},Vec4{0,0,0,0}});
}();

HermiteSpline2D::Mat4x4 HermiteSpline2D::MextR = [](){
    // [0,0,0,0],
    // [0,0,0,0],
    // [0,0,0,lin],
    // [0,1,0,-lin]
    return (Mat4x4 {Vec4{0,0,0,0},Vec4{0,0,0,1},Vec4{0,0,0,0},Vec4{0,0,1,-1}});
}();

HermiteSpline2D::Mat4x4 HermiteSpline2D::MextLT = HermiteSpline2D::MextL.t();
HermiteSpline2D::Mat4x4 HermiteSpline2D::MextRT = HermiteSpline2D::MextR.t();

