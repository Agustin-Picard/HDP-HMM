package utils

import tomasvolker.numeriko.core.functions.average
import tomasvolker.numeriko.core.functions.matMul
import tomasvolker.numeriko.core.functions.transpose
import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D
import tomasvolker.numeriko.core.interfaces.array2d.double.DoubleArray2D
import tomasvolker.numeriko.core.interfaces.array2d.double.elementWise
import tomasvolker.numeriko.core.interfaces.factory.toDoubleArray1D
import tomasvolker.numeriko.core.operations.concatenate
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

fun DoubleArray1D.autocorrelationTime(start: Int = 0, stop: Int = size) =
        this[0 until start].average().let {
            listOf(it).toDoubleArray1D()
                .concatenate(this[start until stop])
                .autocorrelationTime()
        }

fun DoubleArray1D.autocorrelationTime(): Double =
        1 + 2 * sumDouble(1 until size) { i -> autocorrelation(i) } / autocorrelation(0)

fun DoubleArray1D.autocorrelation(delta: Int): Double = mean().let { mean ->
    sumDouble(0 until size - delta) { t ->
        (this[t] - mean) * (this[t + delta] - mean)
    } / (size - delta)
}

fun DoubleArray1D.mean(): Double =
        average()