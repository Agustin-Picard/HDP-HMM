package utils

import hdphmm.mean
import tomasvolker.numeriko.core.functions.average
import tomasvolker.numeriko.core.functions.matMul
import tomasvolker.numeriko.core.functions.transpose
import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D
import tomasvolker.numeriko.core.interfaces.array2d.double.DoubleArray2D
import tomasvolker.numeriko.core.interfaces.array2d.double.elementWise
import tomasvolker.numeriko.core.primitives.squared
import tomasvolker.numeriko.core.primitives.sumDouble
import kotlin.math.pow

fun DoubleArray1D.estimateKurtosis() =
        mean().let { mean ->
            sumDouble(0 until size) { i -> (this[i] - mean).pow(4) } /
                    (sumDouble(0 until size) { i -> (this[i] - mean).squared() }.squared() / size) - 3.0
        }

fun DoubleArray2D.pow(x: Double) =
        singularValueDecomposition().let {
            it.u matMul it.s.elementWise { it.pow(x) } matMul it.v.transpose()
        }

fun DoubleArray1D.unbiased() = this - average()
fun sumVector(indices: IntProgression, selector: (Int) -> DoubleArray1D) =
        indices.asSequence()
                .map { selector(it) }
                .reduce { acc, next -> acc + next }