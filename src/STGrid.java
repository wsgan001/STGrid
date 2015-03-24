import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.TimeUnit;


public class STGrid {
	private static final double R = 6371.0 * 1000;		//the earth's radius (in meters)
	private ArrayList<STPoint> points;
	private double[] spaceBottomLeft;
	private double[] spaceTopRight;
	private long timeBottom;
	private long timeTop;
	private Cell[][][] cells;
	
	
	private double spaceDist;
	private long timeDist;
	private int minPoints;
	
	private List<Cluster> clusters;
	
	
	public static void main (String[] args) {
		try {
			long start = System.currentTimeMillis();
			STGrid stgrid = new STGrid("/Users/jcdmb/Documents/workspace/Project/Julio/STDBSCAN/SCF-CHICAGO-CLEANED.csv", 40, 28, 2);
			stgrid.run();
			long end = System.currentTimeMillis();
			System.out.println("STGrid finished in " + (end - start) + " ms");
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
		
	}
	
	/**
	 * Constructor for a new ST-Grid. It only loads the data from the specified file. The actual algorithm is contained in the method run()
	 * @param csvPath Path to the csv file containing the data to be clustered
	 * @param spaceDist the spatial distance for the generation of the grid (in meters)
	 * @param days the temporal distance for the generation of the grid (in days)
	 * @throws IOException If an error occurred while loading the data from the specified file
	 */
	public STGrid(String csvPath, double spaceDist, long days, int minPoints) throws IOException {
		System.out.println("Reading data");
		long start = System.currentTimeMillis();
		points = CSVLoader.stInitializePoints(csvPath);
		long end = System.currentTimeMillis();
		System.out.println("Finished reading data in " + (end - start) + " ms");
		this.spaceDist = spaceDist;
		this.timeDist = TimeUnit.DAYS.toMillis(days);
		this.minPoints = minPoints;
	}
	
	/**
	 * Constructs a 3-dimensional grid (latitude, longitude and time) big enough to contain the space coordinates and the creation date of every report.
	 * The partition of the grid is given by the parameters specified in the constructor (spaceDist and timeDist). Once the grid is done, every point is 
	 * assigned to the cell whose space and time span contains all three coordinates of the point's report. Finally, all neighboring cells with >= points (each)
	 * are merged to form the clusters.
	 */
	public void run() {
		long beginTime = System.currentTimeMillis();
		constructGrid();
		assignPointsToGrid();
		findClusters();
		ArrayList<Integer> labels = new ArrayList<Integer>();
		for (STPoint point : points) {
			if (point.clusterLabel > 1) {
				labels.add(new Integer(point.clusterLabel));
				Collections.sort(labels);
			}
		}
		for (Integer point : labels) {
			System.out.println("ClusterLabel: "+point);
		}
		long time = System.currentTimeMillis() - beginTime;		
		saveResults("STGRID"+"_"+minPoints+"_"+spaceDist+"_"+timeDist+".csv", 
				"STDBSCAN-Clusters"+"_"+minPoints+"_"+spaceDist+"_"+timeDist+".csv");
		String eval = runAnalysis("STDBSCAN-Analysis"+"_"+minPoints+"_"+spaceDist+"_"+timeDist+".txt", time, spaceDist, timeDist, minPoints, 0, points, clusters.size());
		System.out.println(eval);
		double purity = STGrid.purity(points);
		System.out.println("Purity: "+purity);
		double compression = compression(points);
		System.out.println("Compression: "+compression);
		double entropy = entropyDuplicate(points);
		System.out.println("Entropy: "+entropy);
	}
	
	public List<Cluster> getClusters() {
		return clusters;
	}
	
	public Cell[][][] getGrid() {
		return cells;
	}

	//Constructs the grid by expanding it in every direction on a factor given by the dimensions provided by getGridDimensions().
	private void constructGrid() {
		System.out.println("Constructing grid...");
		constructBoundingBox();
		long start = System.currentTimeMillis();
		int[] dimensions = getGridDimensions();
		int latCells = dimensions[0];
		int lngCells = dimensions[1];
		int timeCells = dimensions[2];
		cells = new Cell[latCells][lngCells][timeCells];
		long end = System.currentTimeMillis();
		System.out.println("Finished constructing grid in " + (end - start) + " ms");
	}
	
	//Reads every STPoint and calculates the borders of the bounding box containing all of them (3-dimensional box: latitude, longitude and time).
		private void constructBoundingBox() {
			System.out.println("Constructing bounding box...");
			long start = System.currentTimeMillis();
			double lowestLat = Double.POSITIVE_INFINITY;
			double highestLat = Double.NEGATIVE_INFINITY;
			double lowestLng = Double.POSITIVE_INFINITY;
			double highestLng = Double.NEGATIVE_INFINITY;
			long lowestTime = Long.MAX_VALUE;
			long highestTime = Long.MIN_VALUE;
			
			for(STPoint point : points) {
				if(point.lat < lowestLat) {
					lowestLat = point.lat;
				} 
				
				if(point.lat > highestLat) {
					highestLat = point.lat;
				}
				
				if(point.lon < lowestLng) {
					lowestLng = point.lon;
				}
				if(point.lon > highestLng) {
					highestLng = point.lon;
				}
				
				if(point.creationDate.getTime() < lowestTime) {
					lowestTime = point.creationDate.getTime();
				} 
				if(point.creationDate.getTime() > highestTime) {
					highestTime = point.creationDate.getTime();
				}
			}
			if(lowestLat == Double.POSITIVE_INFINITY 
					|| highestLat == Double.NEGATIVE_INFINITY 
					|| lowestLng == Double.POSITIVE_INFINITY 
					|| highestLng == Double.NEGATIVE_INFINITY
					|| lowestTime == Long.MAX_VALUE
					|| highestTime == Long.MIN_VALUE) {
				throw new IllegalStateException("Bounding box could not be correctly initialized");
			}
			
			spaceBottomLeft = new double[] {lowestLat, lowestLng};
			spaceTopRight = new double[] {highestLat, highestLng};
			
			//to avoid inaccuracy errors, extend the bounding box by some meters in opposites directions
			//(i.e. pull bottomLeft to the left and further down, and topRight to the right and further up)
			spaceBottomLeft = shiftCoords(spaceBottomLeft, 20, 180);
			spaceBottomLeft = shiftCoords(spaceBottomLeft, 20, -90);
			
			spaceTopRight = shiftCoords(spaceTopRight, 20, 0);
			spaceTopRight = shiftCoords(spaceTopRight, 20, 90);
			
			timeBottom = lowestTime;
			timeTop = highestTime;
			long end = System.currentTimeMillis();
			System.out.println("Finished constructing bounding box in " + (end - start) + " ms.");
		}
	
	//Pre-calculates the dimensions of the grid. Otherwise the grid would have to be initialized using dynamic lists (containing other lists) 
	//instead of multidimensional arrays. 
	private int[] getGridDimensions() {
		System.out.println("Calculating grid dimensions...");
		long start = System.currentTimeMillis();
		double latDistance = getSpaceDistance(spaceBottomLeft, new double[]{spaceTopRight[0], spaceBottomLeft[1]});
		int latCells = (int) (latDistance / spaceDist);
		double endLat = shiftCoords(spaceBottomLeft, (latCells) * spaceDist, 0)[0];
		while(endLat <= spaceTopRight[0]) {
			latCells++;
			endLat = shiftCoords(spaceBottomLeft, (latCells) * spaceDist, 0)[0];
		}
		
		
		double lngDistance = getSpaceDistance(spaceBottomLeft, new double[]{spaceBottomLeft[0], spaceTopRight[1]});
		int lngCells = (int) (lngDistance / spaceDist);
		double endLng = shiftCoords(spaceBottomLeft, (lngCells) * spaceDist, 90)[1];
		while(endLng <= spaceTopRight[1]) {
			lngCells++;
			endLng = shiftCoords(spaceBottomLeft, (lngCells) * spaceDist, 90)[1];
		}
		
		int timeCells = (int) ((timeTop - timeBottom) / timeDist);
		long endTime = timeBottom + (timeCells * timeDist);
		while(endTime <= timeTop) {
			timeCells++;
			endTime = timeBottom + (timeCells * timeDist);
		}
		long end = System.currentTimeMillis();
		System.out.println("Finished calculating grid dimensions in " + (end - start) + " ms");
		System.out.println("Dimensions are: " + latCells + ", " + lngCells + ", " + timeCells);
		
		return new int[]{latCells, lngCells, timeCells};
	}
	
	//Iterates through every STPoint and assigns it to the cell whose spatial and temporal span contains the the spatial and temporal coordinates 
	//of the STPoint being handled.
	//Throws IllegalStateException in case that the grid was not correctly built in one dimension (or more) to contain all STPoints being handled.
	private void assignPointsToGrid() {
		System.out.println("Assigning points to grid...");
		long start = System.currentTimeMillis();
		for(STPoint p : points) {
			int[] coords = getGridCoordinates(p);
			Cell c = cells[coords[0]][coords[1]][coords[2]];
			if(c == null) {
				double startLat = shiftCoords(spaceBottomLeft, coords[0] * spaceDist, 0)[0];
				double endLat = shiftCoords(spaceBottomLeft, (coords[0]+1) * spaceDist, 0)[0];
				double startLng = shiftCoords(spaceBottomLeft, coords[1] * spaceDist, 90)[1];
				double endLng = shiftCoords(spaceBottomLeft, (coords[1]+1) * spaceDist, 90)[1];
				long startTime = timeBottom + (coords[2] * timeDist);
				long endTime = timeBottom + ((coords[2]+1) * timeDist);
				c = new Cell(startLat, endLat, startLng, endLng, startTime, endTime);
				cells[coords[0]][coords[1]][coords[2]] = c;
			}
				c.addPoint(p);
		}
		
		long end = System.currentTimeMillis();
		System.out.println("Finished assigning points to grid in " + (end - start) + " ms");
	}
	
	private int[] getGridCoordinates(STPoint p) {
		int i, j, k;
		for(i = 0; i < cells.length; i++) {
			double startLat = shiftCoords(spaceBottomLeft, i * spaceDist, 0)[0];
			double endLat = shiftCoords(spaceBottomLeft, (i+1) * spaceDist, 0)[0];
			if(startLat <= p.lat && p.lat <= endLat)
				break;
		}
		
		if(i == cells.length)
			throw new IllegalStateException("Grid was not correctly built: a report has a latitude value outside the grid: " + p.lat);
		
		Cell[][] plane = cells[i];
		
		for(j = 0; j < plane.length; j++) {
			double startLng = shiftCoords(spaceBottomLeft, j * spaceDist, 90)[1];
			double endLng = shiftCoords(spaceBottomLeft, (j+1) * spaceDist, 90)[1];
			if(startLng <= p.lon && p.lon <= endLng)
				break;
		}
		
		if(j == plane.length) {
			throw new IllegalStateException("Grid was not correctly built: a report has a longitude value outside the grid: " + p.lon);
		}
		
		Cell[] row = plane[j];
		
		for(k = 0; k < row.length; k++) {
			long startTime = timeBottom + (k * timeDist);
			long endTime = timeBottom + ((k+1) * timeDist);
			if(startTime <= p.creationDate.getTime() && p.creationDate.getTime() <= endTime)
				break;
		}
		
		if(k == row.length) 
			throw new IllegalStateException("Grid was not correctly built: a report has a creation date outsisde the grid.");
		
		return new int[] {i, j, k};
	}
	
	//Iterates through all cells in the grid and for every cell with >= minPoints, a new cluster is grown 
	//(unless said cell already was assigned to another cluster)
	private void findClusters() {
		System.out.println("Clustering points...");
		long start = System.currentTimeMillis();
		clusters = new LinkedList<Cluster>();
		int clusterID = 1;
		for(int i = 0; i < cells.length; i++) {
			for(int j = 0; j < cells[i].length; j++) {
				for(int k = 0; k < cells[i][j].length; k++) {
					Cell c = cells[i][j][k];
					if(c != null && !c.isClustered()) {
						List<STPoint> clusteredPoints = growCluster(i,j,k);
						if(clusteredPoints != null) {
							Cluster cl = new Cluster(clusteredPoints, clusterID++);
							clusters.add(cl);
						}
					}
				}
			}
		}
		long end = System.currentTimeMillis();
		System.out.println("Finished clustering in " + (end - start) + " ms");
	}
	
	//Attempts to grow a new cluster from the cell with the specified array coordinates. 
	//If said cell has not been clustered before and has the sufficient amount of points assigned to it, a new cluster can be expanded.
	//This is done by adding all the points in the current cell to a list and then calling the function recursively for every neighboring cell. 
	private List<STPoint> growCluster(int i, int j, int k) {
		if(i < 0 || j < 0 || k < 0 || i >= cells.length || j >= cells[i].length || k >= cells[i][j].length)
			return null;
			
		Cell c = cells[i][j][k];
		
		if(c == null || c.isClustered() || c.getPointsCount() < minPoints)
			return null;
		
		List<STPoint> clusteredPoints = c.getPoints();
		c.setClustered(true);
		
		List<STPoint> morePoints;
		morePoints = growCluster(i+1, j, k);
		if(morePoints != null)
			clusteredPoints.addAll(morePoints);
		
		morePoints = growCluster(i, j + 1, k);
		if(morePoints != null)
			clusteredPoints.addAll(morePoints);
		
		morePoints = growCluster(i, j, k + 1);
		if(morePoints != null)
			clusteredPoints.addAll(morePoints);
		
		morePoints = growCluster(i-1, j, k);
		if(morePoints != null)
			clusteredPoints.addAll(morePoints);
		
		morePoints = growCluster(i, j - 1, k);
		if(morePoints != null)
			clusteredPoints.addAll(morePoints);
		
		morePoints = growCluster(i, j, k - 1);
		if(morePoints != null)
			clusteredPoints.addAll(morePoints);
		
		return clusteredPoints;
	}
	
	
	
	//Returns an array {latitude, longitude} with the new coordinates in degrees, with a specific distance (in meters) from the origin.
	//The bearing specifies in which direction the origin should be 'shifted'.
	//See http://www.movable-type.co.uk/scripts/latlong.html, function rhumbDestinationPoint()
	//Use 0 as bearing to increase only the latitude, 90 to increase only the longitude
	private double[] shiftCoords(double[] origin, double distance, double bearing) {
		double delta = distance / R;
		double phi1 = Math.toRadians(origin[0]);
		double lambda1 = Math.toRadians(origin[1]);
		double theta = Math.toRadians(bearing);
		
		double dPhi = delta * Math.cos(theta);
		
		double phi2 = phi1 + dPhi;
		
		if(Math.abs(phi2) > Math.PI / 2)
			phi2 = (phi2 > 0) ? Math.PI - phi2 : -Math.PI - phi2;
		
		double dPsi = Math.log(Math.tan(phi2 / 2 + Math.PI / 4) / Math.tan(phi1 / 2 + Math.PI / 4));
		double q = Math.abs(dPsi) > Math.pow(10, -12) ? dPhi / dPsi : Math.cos(phi1);
		
		double dLambda = delta * Math.sin(theta) / q;
		
		double lambda2 = lambda1 + dLambda;
		
		lambda2 = (lambda2 + 3*Math.PI) % (2*Math.PI) - Math.PI;
		
		double[] result = {Math.toDegrees(phi2), Math.toDegrees(lambda2)};
		return result;
	}
	
	private double getSpaceDistance(double[] p1, double[] p2) {
		
		double latn1 = p1[0];
		double longn1 = p1[1];
		double latn2 = p2[0];
		double longn2 = p2[1];
		
		double dLat = Math.toRadians(latn1 - latn2);
		double dLon = Math.toRadians(longn1 - longn2);
		double lat1 = Math.toRadians(latn1);
		double lat2 = Math.toRadians(latn2);
		
		double a = Math.sin(dLat / 2) * Math.sin(dLat / 2) + Math.sin(dLon / 2) * Math.sin(dLon / 2) * Math.cos(lat1) * Math.cos(lat2);
		a = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
		a = R * a;
		
		return a;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Test methods (help to find the problem in case some IllegalStateException is thrown)
	
	
	//Method for testing the generation of the grid
	@SuppressWarnings("unused")
	private void testGridArray() {
		//Test latitudes:
		for(int i = 0; i < cells.length; i++) {
			Cell[][] plain = cells[i];
			for(int j = 0; j < plain.length; j++) {
				Cell[] row = plain[j];
				for(int k = 0; k < row.length - 1; k++) {
					Cell c1 = row[k];
					Cell c2 = row[k+1];
					if(Math.abs(c1.getStartLat() - c2.getStartLat()) > 0.001)
						throw new IllegalStateException("Start latitudes in the same xy-row do not match");
					if(Math.abs(c1.getEndLat() - c2.getEndLat()) > 0.001)
						throw new IllegalStateException("End latitudes in the same xy-row do not match");
					if(Math.abs(c1.getStartLng() - c2.getStartLng()) > 0.001)
						throw new IllegalStateException("Start longitudes in the same xy-row do not match");
					if(Math.abs(c1.getEndLng() - c2.getEndLng()) > 0.001)
						throw new IllegalStateException("End longitudes in the same xy-row do not match");
				}
			}
		}
		
		for(int i = 0; i < cells.length; i++) {
			for(int j = 0; j < cells[i].length - 1; j++) {
				for(int k = 0; k < cells[i][j].length; k++) {
					Cell c1 = cells[i][j][k];
					Cell c2 = cells[i][j + 1][k];
					if(Math.abs(c1.getStartLat() - c2.getStartLat()) > 0.001)
						throw new IllegalStateException("Start latitudes in the same xz-row do not match");
					if(Math.abs(c1.getEndLat() - c2.getEndLat()) > 0.001)
						throw new IllegalStateException("End latitudes in the same xz-row do not match");
					if(Math.abs(c1.getStartTime() - c2.getStartTime()) > 0.001)
						throw new IllegalStateException("Start times in the same xz-row do not match");
					if(Math.abs(c1.getEndTime() - c2.getEndTime()) > 0.001)
						throw new IllegalStateException("End times in the same xz-row do not match");
				}
			}
		}
		
		for(int i = 0; i < cells.length - 1; i++) {
			for(int j = 0; j < cells[i].length; j++) {
				for(int k = 0; k < cells[i][j].length; k++) {
					Cell c1 = cells[i][j][k];
					Cell c2 = cells[i + 1][j][k];
					if(Math.abs(c1.getStartLng() - c2.getStartLng()) > 0.001)
						throw new IllegalStateException("Start longitudes in the same yz-row do not match");
					if(Math.abs(c1.getEndLng() - c2.getEndLng()) > 0.001)
						throw new IllegalStateException("End longitudes in the same yz-row do not match");
					if(Math.abs(c1.getStartTime() - c2.getStartTime()) > 0.001)
						throw new IllegalStateException("Start times in the same yz-row do not match");
					if(Math.abs(c1.getEndTime() - c2.getEndTime()) > 0.001)
						throw new IllegalStateException("End times in the same yz-row do not match");
				}
			}
		}
	}
	
	//Test consistency of cell dimensions
	@SuppressWarnings("unused")
	private void testCells() {
		int nonEmptyCells = 0;
		int points = 0;
		int maxPointsProCell = 0;
		for(int i = 0; i < cells.length; i++) {
			Cell[][] plain = cells[i];
			for(int j = 0; j < plain.length; j++) {
				Cell[] row = plain[j];
				for(int k = 0; k < row.length; k++) {
					Cell c = row[k];
					if(c != null ) {
						System.out.println("Cell (" + i + ", " + j + ", " + k + ") contains: " + c.getPointsCount() + " points");
						nonEmptyCells++;
						points += c.getPointsCount();
						if(c.getPointsCount() > maxPointsProCell)
							maxPointsProCell = c.getPointsCount();
					}
				}
			}
		}
		System.out.println("Total non empty cells: " + nonEmptyCells + ", total reports: " + points + ", max points in a cell: " + maxPointsProCell);
	}
	
	public void printBoundingBox() {
		System.out.println("-----------------------------------");
		System.out.println("Bounding box:");
		
		double[] spaceBottomRight = new double[]{spaceBottomLeft[0], spaceTopRight[1]};
		double[] spaceTopLeft = new double[]{spaceTopRight[0], spaceBottomLeft[1]};
		
		System.out.println("\tspaceBottomLeft: " + spaceBottomLeft[0] + ", " + spaceBottomLeft[1]);
		System.out.println("\tspaceBottomRight: " + spaceBottomRight[0] + ", " + spaceBottomRight[1]);
		System.out.println("\tspaceTopLeft: " + spaceTopLeft[0] + ", " + spaceTopLeft[1]);
		System.out.println("\tspaceTopRight: " + spaceTopRight[0] + ", " + spaceTopRight[1]);
		System.out.println("\ttimeBottom: " + timeBottom + ", timeTop: " + timeTop);
		
		double latSpan = getSpaceDistance(spaceBottomLeft, spaceTopLeft);
		System.out.println("\t\tLatitude span: " + latSpan + " m");
		System.out.println("\t\t#lat cells: " + (latSpan / spaceDist));
		
		double lngSpan = getSpaceDistance(spaceBottomLeft, spaceBottomRight);
		System.out.println("\t\tLongitude span: " + lngSpan + " m");
		System.out.println("\t\t#lng cells: " + (lngSpan / spaceDist));
		
		double timeSpan = Math.abs(timeTop - timeBottom);
		System.out.println("\t\tTime span: " + timeDist + " ms");
		System.out.println("\t\t#time cells: " + (double) (timeSpan / timeDist));
		System.out.println("-----------------------------------");
				

	}
	
	
	//Method for testing the generated clusters.
	@SuppressWarnings("unused")
	private void testClusters() {
		for(Cluster cl : clusters) {
			System.out.println("-----------------------------------");
			System.out.println("Cluster " + cl.clusterLabel + ": ");
			List<STPoint> clPoints = cl.points;
			for(int i = 0; i < clPoints.size(); i++) {
				STPoint currentPoint = clPoints.get(i);
				for(int j = i+1; j < clPoints.size(); j++) {
					STPoint nextPoint = clPoints.get(j);
					double spaceDistance = getSpaceDistance(new double[]{currentPoint.lat, currentPoint.lon}, new double[]{nextPoint.lat, nextPoint.lon});
					double timeDistance = Math.abs(currentPoint.creationDate.getTime() - nextPoint.creationDate.getTime());
					if(spaceDistance > spaceDist || timeDistance > timeDist) {
						System.err.println("p" + currentPoint.index + " - p" + nextPoint.index + ": spaceDist: " + spaceDistance + ", timeDist: " + timeDistance);
					} else {
						System.out.println("p" + currentPoint.index + " - p" + nextPoint.index + ": spaceDist: " + spaceDistance + ", timeDist: " + timeDistance);
					}
				}
			}
		}
	}
	
	/**
	 * Save clustering Results to two CSV files.
	 * 1) pointsPath = Original Data with extra column containing cluster label
	 * 2) clustersPath = Cluster center coordinates and cluster size
	 */
	protected void saveResults(String pointsPath, String clustersPath) {
		try {
			PrintWriter writer = new PrintWriter(clustersPath, "UTF-8");
			writer.write(Cluster.CSV_Header());
			for (Cluster cluster : clusters) {
				if (cluster.clusterLabel > 0)
					writer.write(cluster.toCSV());
			}
			writer.close();

			writer = new PrintWriter(pointsPath, "UTF-8");
			writer.write(STPoint.CSV_Header());
			for (Cluster cluster : clusters) {
				for (STPoint point : cluster.points) {
					writer.write(point.toCSV());
				}
			}
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static String runAnalysis(String resultsPath, long time, double eps, long eps2, int minPointThreshHold, double deltaE, List<STPoint> points, int clusterLabel) {
		 double[] metrics = evalMetrics(points);
		 double precision = metrics[0];
		 double recall = metrics[1];
		 double fmeasure = metrics[2];
		 double total = metrics[3];
		 double totalInCluster = metrics[4];
		 double tp = metrics[5];
		 double fp = metrics[6];
		 double tn = metrics[7];
		 double fn = metrics[8];
		 
		 double classErrorTrue = fn/(tp+fn);
		 double classErrorFalse = fp/(tn+fp);
		
		String analysis = "\nThe Algorithm took " + time
					+ " milliseconds to complete for "+points.size()+" data points"+"\n";
		analysis += "#Clusters found: " + clusterLabel+"\n";
		analysis += (deltaE > 0)?"Using Leventhstein Similarity\n":"";
		analysis += "Configuration - eps1: "+eps+"m, eps2: "+eps2/86400000+" days, minpts: "+minPointThreshHold+", deltaE: "+deltaE+"\n";
		analysis += "#Points: "+total+". #Points in clusters: "+totalInCluster+"\n";
		analysis += "TP: "+tp+", FP: "+fp+", TN: "+tn+", FN: "+fn+"\n";
		analysis += "Class.Error.True: "+classErrorTrue+", Class.Error.False: "+classErrorFalse+"\n";
		analysis += "Precision: "+precision+". Recall: "+recall+". F-Measure: "+fmeasure+"\n";
		analysis += "Runtime: "+time/1000+" seconds";
		
		if (resultsPath != null) {
			PrintWriter writerCenter;
			try {
				writerCenter = new PrintWriter(resultsPath, "UTF-8");
				writerCenter.println(analysis);
				writerCenter.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return analysis;
	}
	
	public static double[] evalMetrics(List<STPoint> points) {
		int total=0;
		int totalInCluster=0;
		double fp = 0;
		double tp = 0;
		double tn = 0;
		double fn = 0;
		
		for (STPoint p : points) {
			total++;
			if (p.isDuplicate && p.clusterLabel>0) {tp++; totalInCluster++;}
			else if (!p.isDuplicate && p.clusterLabel>0) {fp++;totalInCluster++;}
			else if (p.isDuplicate && !(p.clusterLabel>0)) fn++;
			else if (!p.isDuplicate && !(p.clusterLabel>0)) tn++;
		}
		
		double precision = tp/(tp + fp);
		double recall = tp/(tp+fn);
		double fmeasure = (2*precision*recall)/(precision+recall);
		
		return new double[] {precision,recall,fmeasure, total, totalInCluster,tp,fp,tn,fn};
	}
	
	public static int clusterPurity(Cluster cluster) {
		int trueLabel=0;
		int falseLabel=0;
		for (STPoint point : cluster.points) {
			if (point.isDuplicate)
				trueLabel++;
			else
				falseLabel++;
		}
		return Math.max(trueLabel, falseLabel);
	}
	
	public static <T extends Comparable<? super T>> List<T> asSortedList(
			Collection<T> c) {
		List<T> list = new ArrayList<T>(c);
		java.util.Collections.sort(list);
		return list;
	}
	
	@SuppressWarnings("unchecked")
	private static void addPoint2Clustermap(
			@SuppressWarnings("rawtypes") HashMap multiMap, STPoint p,
			int clusterLabel) {
		List<STPoint> list;
		if (multiMap.containsKey(p.clusterLabel)) {
			list = (List<STPoint>) multiMap.get(p.clusterLabel);
			list.add(p);
		} else {
			list = new ArrayList<STPoint>();
			list.add(p);
			multiMap.put(p.clusterLabel, list);
		}
	}
	public static ArrayList<Cluster> generateClusters(ArrayList<STPoint> points) {
		HashMap<Integer, List<STPoint>> clusterTable = new HashMap<Integer, List<STPoint>>();
		ArrayList<Cluster> clusters = new ArrayList<Cluster>();
		for (STPoint p : points) {
			addPoint2Clustermap(clusterTable, p, p.clusterLabel);
		}

		for (Integer key : asSortedList(clusterTable.keySet())) {
			List<STPoint> list = clusterTable.get(key);
			Cluster cluster = new Cluster(list, key);
			clusters.add(cluster);
		}
		return clusters;
	}
	
	public static double purity(ArrayList<STPoint> points) {
		//private static ArrayList<Cluster> generateClusters(ArrayList<STPoint> points) {
		ArrayList<Cluster> clusters = generateClusters(points);
		int purityInt = 0;
		int size = 0;
		for (Cluster cluster : clusters) {
			purityInt = (int) (purityInt + clusterPurity(cluster));
			size = size + cluster.points.size();
		}
		return (double)((double)purityInt/(double)size);
	}
	
	public static double compression(ArrayList<STPoint> points) {
		int originalSize = points.size();
		int clustersAmount = 0;
		int outliersAmount = 0;
		for (STPoint point : points) {
			if (point.clusterLabel>clustersAmount)
				clustersAmount = point.clusterLabel;
			if (point.clusterLabel==-1)
				outliersAmount = outliersAmount+1;
		}
		return (double)((double)(outliersAmount+clustersAmount)/(double)originalSize);
	}
	
	public static double entropyDuplicate(ArrayList<STPoint> points) {
		ArrayList<Boolean> labels = new ArrayList<Boolean>();
		double duplicateCount = 0;
		for (STPoint point : points) {
			labels.add(point.isDuplicate);
			labels.add(point.clusterLabel>0);
			if (point.isDuplicate)
				duplicateCount++;
				
		}
		double nonduplicateCount = points.size() - duplicateCount;
		double total = points.size();
		
		double p1 = (nonduplicateCount/total);
		double l1 = Math.log(nonduplicateCount/total);
		double p2 = (duplicateCount/total);
		double l2 = Math.log(duplicateCount/total);
		
		return -(p1*l1)-(p2*l2);
	}
}
