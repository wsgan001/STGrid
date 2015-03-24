
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.concurrent.TimeUnit;


/**
 * Licensed to Creative Commons ShareAlike 4.0 International (CC BY-SA 4.0)
 * */

/**
 * The class models the spatio-temporal points to be clustered.
 * It possses spatio, temporal and semantic parameters as well as
 * parameters related to the clustering itself. Methods for comparing
 * points w.r.t to these parameters are included in this class.
 * 
 * @author De Melo Borges
 */
public class STPoint {
	public int index;
	public boolean isDuplicate = false;
	public boolean visited = false;
	// SPATIAL ATTRIBUTES:
	public double lat;
	public double lon;
	public double x;
	public double y;
	public double z;
	// TEMPORAL ATTRIBUTE
	public Date creationDate;
	//SEMANTIC ATTRIBUTES
	public String text;
	public String url;
	// CLUSTERING ATTRIBUTES
	public boolean isNoise = false;
	public int clusterLabel = -1;
	// Radius of the earth in meters
	public static final double R = 6371.0 * 1000;
	//SEMANTIC ATTRIBUTE
	public STPoint(double lat, double lon, boolean isDuplicate,
			Date creationDate, int index, String text, String url) {
		// super();
		this.index = index;
		this.lat = lat;
		this.lon = lon;
		this.isDuplicate = isDuplicate;
		this.creationDate = creationDate;
		this.text = text;
		this.url = url;
		gisToEuclidian(lat, lon);
	}
	
	public STPoint(double lat, double lon) {
		this.lat = lat;
		this.lon = lon;
	}

	/**
	 * Get the difference between two dates in days
	 * 
	 * @param t - the comparing point.
	 * @return the difference value in days.
	 */
	public long temporalDistanceFromPoint(STPoint t) {
		Date date1 = creationDate;
		Date date2 = t.creationDate;
		TimeUnit timeUnit = TimeUnit.DAYS;
		long diffInMillies = date2.getTime() - date1.getTime();
		return Math.abs(timeUnit.convert(diffInMillies, TimeUnit.MILLISECONDS));
	}

	
	/**
	 * Calls the suitable method for calculating the spatial
	 * distance between two points. It depends on STDBSCAN.isEuclidian.
	 * If this is flagged to true, the KD-Tree of the class will 
	 * calculate distances based on the euclidian paremeters, so the 
	 * method should delivers the euclidian based distance as well.
	 * Otherwise delivers the haversine distance based on the coordinates.  
	 * 
	 * @param t - the comparing point.
	 * @return the spatial distance in meters.
	 */
	public double spatialDistanceFromPoint(STPoint t) {
			return haversineDist(t);
	}

	/**
	 * Euclidian Distance
	 * 
	 * @param t - the comparing point.
	 * @return the spatial distance in meters.
	 */
	public double euclidianDistance(STPoint t) {
		return Math.sqrt(Math.pow(x - t.x, 2) + Math.pow(y - t.y, 2)
				+ Math.pow(z - t.z, 2));
	}

	/**
	 * Haversine Distance based on coordinates.
	 * 
	 * @param t - the comparing point.
	 * @return the spatial distance in meters.
	 */
	public double haversineDist(STPoint t) {
		double lat1 = this.lat;
		double lng1 = this.lon;
		double lat2 = t.lat;
		double lng2 = t.lon;
		double dLat = Math.toRadians(lat2 - lat1);
		double dLng = Math.toRadians(lng2 - lng1);
		double sindLat = Math.sin(dLat / 2);
		double sindLng = Math.sin(dLng / 2);
		double a = Math.pow(sindLat, 2) + Math.pow(sindLng, 2)
				* Math.cos(Math.toRadians(lat1))
				* Math.cos(Math.toRadians(lat2));
		double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
		double dist = (R * c);
		return dist;
	}

	/**
	 * Utility method to transform the geographical 
	 * coordinates to euclidian coordinates.
	 * 
	 * @param lat
	 * @param lon
	 */
	private void gisToEuclidian(double lat, double lon) {
		lat *= Math.PI / 180.0;
		lon *= Math.PI / 180.0;

		this.x = R * Math.cos(lat) * Math.cos(lon);
		this.y = R * Math.cos(lat) * Math.sin(lon);
		this.z = R * Math.sin(lat);
	}
	
//	public double levenshteinSim(STPoint t) {
//		//1 - d(str1,str2) / max(A,B)
//		double levDist = (double) StringUtils.getLevenshteinDistance(this.text, t.text);
//		return 1 - (levDist)/Math.max(this.text.length(), t.text.length());
//	}

	
	/**
	 * Euclidian coordinates as array.
	 * Used as key entry by the KD-Tree
	 * @return array containing euclidian coordinates.
	 */
	public double[] toDouble() {
		double[] xyz = { x, y, z };
		return xyz;
	}

	/**
	 * Utility method needed for exporting points to CSV
	 * @return the CSV String, separated by commas
	 */
	public String toCSV() {
		return index+","+dateFormated()+","+lat+","+lon+","+x+","+y+","+z+","+url+","+isDuplicate+","+clusterLabel+"\n";
		
	}
	
	/**
	 * Utility method needed to delivers the date (creation date)
	 * in suitable format for file writing
	 * @return formatted date as String
	 */
	private String dateFormated() {
		SimpleDateFormat simpleDateFormat = new SimpleDateFormat(
				"yyyy-MM-dd HH:mm:ss");
		return simpleDateFormat.format(this.creationDate);
		
	}
	
	/**
	 * The Header for CSV File 
	 * @return
	 */
	public static String CSV_Header() {
		return "index,creationDate,lat,lon,x,y,z,url,duplicate,clusterLabel\n";
	}
	
	public boolean equals(STPoint t) {
		return (lat == t.lat) && (lon == t.lon);
	}

	public String toString() {
		return "Created at: " + this.creationDate.toString()
				+ ". Coordinates: " + this.lat + ", " + this.lon;
	}
}
