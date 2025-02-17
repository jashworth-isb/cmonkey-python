package controllers

import play.api.libs.json._
import java.io._
import java.util.regex._
import scala.collection.mutable.HashMap

case class Annotation(motifNum: Int, reverse: Boolean, position: Int, gene: String, pvalue: Double)
case class GeneAnnotations(gene: String, annotations: Seq[Annotation])
case class MotifInfo(motifNum: Int, evalue: Double, pssm: Array[Array[Float]], annotations: Array[Annotation])
case class Snapshot(rows: Map[Int, List[String]], columns: Map[Int, List[String]],
                    residuals: Map[Int, Double],
                    motifs: Map[String, Map[Int, Array[MotifInfo]]]) {
}

object SnapshotReader {
  val BaseResultsFileName = "%d-results.json"
  val JsonFilePattern = Pattern.compile("(\\d+)-results.json")
}
class SnapshotReader(OutDirectory: File, Synonyms: SynonymsMap) {
  import SnapshotReader._

  implicit object SnapshotFormat extends Format[Snapshot] {

    private def readMotifInfos(stMotifs: Seq[JsValue]) = {
      val motifInfos = new java.util.ArrayList[MotifInfo]
      for (motif <- stMotifs) {
        val motifObj = motif.asInstanceOf[JsObject]
        val pssm = motifObj.value("pssm").as[Array[Array[Float]]]
        val evalue = if (motifObj.keys.contains("evalue")) (motifObj \ "evalue").as[Double]
                     else 0.0
        val motifNum = (motifObj \ "motif_num").as[Int]
        val annotations = if (motifObj.keys.contains("annotations")) {
          val annots = (motifObj \ "annotations").asInstanceOf[JsArray]
          val annotArr = new Array[Annotation](annots.value.length)
          for (i <- 0 until annotArr.length) {
            val current = annots(i)
            annotArr(i) = Annotation(motifNum,
                                     (current \ "reverse").as[Boolean],
                                     (current \ "position").as[Int],
                                     (current \ "gene").as[String],
                                     (current \ "pvalue").as[Double])
          }
          annotArr
        } else Array[Annotation]()
        motifInfos.add(MotifInfo(motifNum, evalue, pssm, annotations))
      }
      motifInfos.toArray(new Array[MotifInfo](0))
    }

    def reads(json: JsValue): Snapshot = {
      val rows = (json \ "rows").as[JsObject]
      val cols = (json \ "columns").as[JsObject]
      val residuals = (json \ "residuals").as[JsObject]
      val motifsVal = (json \ "motifs")
      val clusterRows = new HashMap[Int, List[String]]
      for (field <- rows.fields) {
        clusterRows(field._1.toInt) = field._2.as[List[String]].map(str => Synonyms(str))
      }

      val clusterCols = new HashMap[Int, List[String]]
      for (field <- cols.fields) {
        clusterCols(field._1.toInt) = field._2.as[List[String]]
      }

      val clusterResiduals = new HashMap[Int, Double]
      for (field <- residuals.fields) {
        clusterResiduals(field._1.toInt) = field._2.asInstanceOf[JsNumber].value.doubleValue
      }

      //val clusterMotifs = new HashMap[Int, Map[String, Array[MotifInfo]]]
      val seqTypeMotifs = new HashMap[String, Map[Int, Array[MotifInfo]]]
      try {
        motifsVal match {
          case motifs:JsObject =>
            for (field <- motifs.fields) {
              val seqType = field._1
              val clusterObj = field._2.as[JsObject] // Map[Int, Array[MotifInfo]]
              //val seqTypeObj = field._2.as[JsObject]
              //val seqTypeMotifs = new HashMap[String, Array[MotifInfo]]
              val clusterMotifs = new HashMap[Int, Array[MotifInfo]]

              // iterate over the clusters, which are the keys
              for (cluster <- clusterObj.keys) {
                // an array of motif objects (motif_num, evalue, annotations, pssm)
                // annotations are tuples of (gene, position, pvalue, reverse)
                val stResult = clusterObj \ cluster
                val motifInfos = (stResult \ "motif-info")
                  //seqTypeMotifs(seqType) = if (motifInfos.isInstanceOf[JsArray]) {
                  //  readMotifInfos(motifInfos.asInstanceOf[JsArray].value)
                  clusterMotifs(cluster.toInt) = if (motifInfos.isInstanceOf[JsArray]) {
                    readMotifInfos(motifInfos.asInstanceOf[JsArray].value)
                  } else {
                    Array.ofDim[MotifInfo](0)
                  }
              }
              //clusterMotifs(field._1.toInt) = seqTypeMotifs.toMap
              seqTypeMotifs(seqType) = clusterMotifs.toMap
            }
          case _ =>
            println("no motif values found")
        }
      } catch {
        case e =>
          e.printStackTrace
          println("\nNo motifs found !!!")
      }
      Snapshot(clusterRows.toMap, clusterCols.toMap, clusterResiduals.toMap, seqTypeMotifs.toMap)
    }

    def writes(snapshot: Snapshot): JsValue = JsUndefined("TODO")
  }

  def readSnapshot(iteration: Int) : Option[Snapshot] = {
    val pathname = (OutDirectory + "/" + BaseResultsFileName).format(iteration)
    printf("Reading snapshot: %s\n", pathname)
    val infile = new File(pathname)
    if (infile.exists) {
      val in = new BufferedReader(new FileReader(infile))
      val buffer = new StringBuilder
      var line = in.readLine
      while (line != null) {
        buffer.append(line)
        line = in.readLine
      }
      in.close
      Some(play.api.libs.json.Json.parse(buffer.toString).as[Snapshot])
    } else {
      printf("File '%s' does not exist !\n", infile.getName)
      None
    }
  }
}
