/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Owen Jessop
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include <config.h>

/* Some standard headers. */
#include <math.h>
#include <stdio.h>

/* Local headers. */
#include "kernel_gravity.h"

#define numPoints 1000

/**
 * @brief The softened potential of a unit mass at distance r.
 *
 * @param r The distance.
 * @param H The softening length (kernel support).
 */
static double potential(const float r, const float H) {
  if (r >= H) return -1. / r;
  return kernel_grav_pot_eval(r / H) / H;
}

/**
 * @brief Checks kernel_grav_pot_dH_eval() against a finite-difference
 * derivative of the softened potential with respect to the softening length.
 */
int main(int argc, char *argv[]) {

  const float H = 1.f;
  const float dH = 1e-2f;

  for (int k = 0; k <= numPoints; ++k) {

    const float u = (float)k / numPoints;
    const float r = u * H;

    const double dphi_dH_fd =
        (potential(r, H + dH) - potential(r, H - dH)) / (2. * dH);
    const double dphi_dH = kernel_grav_pot_dH_eval(u) / (H * H);

    if (fabs(dphi_dH - dphi_dH_fd) > 1e-3) {
      printf("u= %f dphi/dH= %e finite difference= %e\n", u, dphi_dH,
             dphi_dH_fd);
      printf("Invalid derivative of the softened potential!\n");
      return 1;
    }
  }

  /* The potential is -1/r outside the kernel for any H */
  if (fabsf(kernel_grav_pot_dH_eval(1.f)) > 1e-5f) {
    printf("dphi/dH= %e at the kernel edge, should be 0!\n",
           kernel_grav_pot_dH_eval(1.f));
    return 1;
  }

  /* At the centre, phi = W_pot(0) / H so dphi/dH = -W_pot(0) / H^2 */
  if (fabsf(kernel_grav_pot_dH_eval(0.f) + kernel_grav_pot_eval(0.f)) > 1e-5f) {
    printf("dphi/dH= %e at the centre, should be %e!\n",
           kernel_grav_pot_dH_eval(0.f), -kernel_grav_pot_eval(0.f));
    return 1;
  }

  printf("\nAll values are consistent\n");
  return 0;
}
