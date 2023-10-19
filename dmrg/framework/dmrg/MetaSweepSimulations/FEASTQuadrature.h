/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef FEAST_QUADRATURE
#define FEAST_QUADRATURE

#include <stdexcept>
#include <vector>

namespace FeastHelper {

class QuadraturePoint {
public:
  /** @brief Class constructor */
  QuadraturePoint(double node, double weight)
    : quadratureNode(node), quadratureWeight(weight) {}
  // Struct members
  double quadratureNode;
  double quadratureWeight;
};

/**
 * @brief Function that generates the quadrature points
 * Inlined to avoid multiple definition problem.
 */
inline std::vector<QuadraturePoint> getQuadraturePoints(int nQuad)
{
  std::vector<QuadraturePoint> ret;
  switch(nQuad) {
    case 2:
      ret = {QuadraturePoint(-0.5773502691896257645092, 1.0),
             QuadraturePoint( 0.5773502691896257645092, 1.0)};
      break;
    case 3:
      ret = {QuadraturePoint(-0.7745966692414833770359, 0.5555555555555555555556),
             QuadraturePoint( 0, 0.8888888888888888888889),
             QuadraturePoint( 0.7745966692414833770359, 0.555555555555555555556)};
      break;
    case 4:
      ret = {QuadraturePoint(-0.861136311594052575224, 0.3478548451374538573731),
             QuadraturePoint(-0.3399810435848562648027, 0.6521451548625461426269),
             QuadraturePoint( 0.3399810435848562648027, 0.6521451548625461426269),
             QuadraturePoint( 0.861136311594052575224, 0.3478548451374538573731)};
      break;
    case 5:
      ret = {QuadraturePoint(-0.9061798459386639927976, 0.2369268850561890875143),
             QuadraturePoint(-0.5384693101056830910363, 0.4786286704993664680413),
             QuadraturePoint( 0, 0.5688888888888888888889),
             QuadraturePoint( 0.5384693101056830910363, 0.4786286704993664680413),
             QuadraturePoint( 0.9061798459386639927976, 0.2369268850561890875143)};
      break;
    case 6:
      ret = {QuadraturePoint(-0.9324695142031520278123, 0.1713244923791703450403),
             QuadraturePoint(-0.661209386466264513661, 0.3607615730481386075698),
             QuadraturePoint(-0.2386191860831969086305, 0.4679139345726910473899),
             QuadraturePoint( 0.238619186083196908631, 0.46791393457269104739),
             QuadraturePoint( 0.661209386466264513661, 0.3607615730481386075698),
             QuadraturePoint( 0.9324695142031520278123, 0.1713244923791703450403)};
      break;
    case 7:
      ret = {QuadraturePoint(-0.9491079123427585245262, 0.1294849661688696932706),
             QuadraturePoint(-0.7415311855993944398639, 0.2797053914892766679015),
             QuadraturePoint(-0.4058451513773971669066, 0.38183005050511894495),
             QuadraturePoint( 0, 0.417959183673469387755),
             QuadraturePoint( 0.4058451513773971669066, 0.38183005050511894495),
             QuadraturePoint( 0.7415311855993944398639, 0.279705391489276667901),
             QuadraturePoint( 0.9491079123427585245262, 0.129484966168869693271)};
      break;
    case 8:
      ret = {QuadraturePoint(-0.9602898564975362316836, 0.1012285362903762591525),
             QuadraturePoint(-0.7966664774136267395916, 0.2223810344533744705444),
             QuadraturePoint(-0.5255324099163289858177, 0.313706645877887287338),
             QuadraturePoint(-0.1834346424956498049395, 0.3626837833783619829652),
             QuadraturePoint( 0.1834346424956498049395, 0.3626837833783619829652),
             QuadraturePoint( 0.5255324099163289858177, 0.313706645877887287338),
             QuadraturePoint( 0.7966664774136267395916, 0.222381034453374470544),
             QuadraturePoint( 0.9602898564975362316836, 0.1012285362903762591525)};
      break;
    case 24:
      ret = { QuadraturePoint(-0.99518721999702136018, 0.0123412297999871995468),
              QuadraturePoint(-0.9747285559713094981984, 0.0285313886289336631813),
              QuadraturePoint(-0.9382745520027327585237, 0.0442774388174198061686),
              QuadraturePoint(-0.8864155270044010342132, 0.05929858491543678074637),
              QuadraturePoint(-0.820001985973902921954, 0.07334648141108030573403),
              QuadraturePoint(-0.7401241915785543642438, 0.0861901615319532759172),
              QuadraturePoint(-0.648093651936975569253, 0.0976186521041138882699),
              QuadraturePoint(-0.545421471388839535658, 0.1074442701159656347826),
              QuadraturePoint(-0.4337935076260451384871, 0.1155056680537256013533),
              QuadraturePoint(-0.3150426796961633743868, 0.1216704729278033912045),
              QuadraturePoint(-0.191118867473616309159, 0.1258374563468282961214),
              QuadraturePoint(-0.064056892862605626085, 0.1279381953467521569741),
              QuadraturePoint(0.064056892862605626085, 0.1279381953467521569741),
              QuadraturePoint(0.1911188674736163091586, 0.1258374563468282961214),
              QuadraturePoint(0.3150426796961633743868, 0.121670472927803391204),
              QuadraturePoint(0.4337935076260451384871, 0.115505668053725601353),
              QuadraturePoint(0.5454214713888395356584, 0.1074442701159656347826),
              QuadraturePoint(0.6480936519369755692525, 0.09761865210411388827),
              QuadraturePoint(0.7401241915785543642438, 0.0861901615319532759172),
              QuadraturePoint(0.820001985973902921954, 0.073346481411080305734),
              QuadraturePoint(0.8864155270044010342132, 0.0592985849154367807464),
              QuadraturePoint(0.9382745520027327585237, 0.044277438817419806169),
              QuadraturePoint(0.9747285559713094981984, 0.0285313886289336631813),
              QuadraturePoint(0.99518721999702136018, 0.0123412297999871995468)};
      break;
    default:
      throw std::runtime_error("Number of quadrature points not supported");
  }
  return ret;
}

} // namespace FeastHelper

#endif // FEAST_QUADRATURE
