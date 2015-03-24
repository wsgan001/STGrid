

import java.util.ArrayList;
import java.util.List;

/**
 * Licensed to Creative Commons ShareAlike 4.0 International (CC BY-SA 4.0)
 * */

/**
 * This class represents a cluster containing the STPoints.
 * Objects of this class are created after the clustering 
 * results of STDBSCAN.java. 
 * 
 * @author Julio De Melo Borges
 */
public class Cluster {

	public ArrayList<STPoint> points = new ArrayList<STPoint>();
	public int clusterLabel;
	
	public Cluster(List<STPoint> list, int clusterLabel) {
		this.points = (ArrayList<STPoint>) list;
		this.clusterLabel = clusterLabel;
		for (STPoint point : list) {
			point.clusterLabel = clusterLabel;
		}
	}

	public int size() {
		return points.size();
	}

	public double latCenter() {
		double latSum = 0;
		for (STPoint point : points) {
			latSum += point.lat;
		}
		return latSum/size();
	}
	
	public double lonCenter() {
		double lonSum = 0;
		for (STPoint point : points) {
			lonSum += point.lon;
		}
		return lonSum/size();
	}

	public static String CSV_Header() {
		return "label,lat,lon,size\n";
	}

	public String toCSV() {
		return clusterLabel+","+latCenter()+","+lonCenter()+","+size()+"\n";
	}
	
}
