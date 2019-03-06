package utils

import java.io.File


data class Trace(val label: Int, val timestamps: List<Double>, val data: List<Double>)

class SegmentedTracesReader(val filename: String) {

    fun readTraces(): List<Trace> {
        val traceList = mutableListOf<Trace>()
        val joinedData = mutableListOf<Pair<List<Double>,Int>>()
        val classList = readClasses().toList()
        val classCount = classList.distinct().size

        readData().forEachIndexed { index, data ->
            joinedData.add(data to classList[index])
        }

        for (i in 0 until classCount) {
            val filteredData = joinedData.filter { it.second == i }.map { it.first }
            traceList.add(Trace(i, filteredData[0], filteredData[1]))
        }

        return traceList
    }

    private fun readData() =
            File(filename).useLines {
                it.map {
                    it.trim()
                    it.split(',').map { it.toDouble() }
                }
            }

    private fun readClasses() =
            File("$filename.seq")
                .useLines {
                    it.map {
                        it.trim()
                            .toInt()
                    }
                }
}