package utils

import org.apache.commons.math3.special.Erf
import org.apache.commons.math3.special.Gamma

fun gammaFunction(x: Double): Double =
        Gamma.gamma(x)

fun logGammaFunction(x: Double) =
        Gamma.logGamma(x)

fun erfinvFunction(x: Double) =
        Erf.erfInv(x)

fun Int.factorial(): Int =
    if (this > 0) this * factorial() else 1
