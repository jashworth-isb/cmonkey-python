package controllers

object Formatter {
  def formatEvalues(motifInfos: Map[String, Map[Int, Array[MotifInfo]]], cluster: Int): String = {
    val comps = motifInfos.keys.map { seqType =>
      val clusterMotifInfos = motifInfos(seqType)
      val evals = if (clusterMotifInfos.contains(cluster)) {
        clusterMotifInfos(cluster).map { motifInfo => 
          "%.2e".format(motifInfo.evalue)
        }.mkString(",")
      } else ""
      "(%s -> [%s])".format(seqType, evals)
    }
    comps.mkString("<br>")
  }
}
