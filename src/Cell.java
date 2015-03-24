import java.util.ArrayList;
import java.util.List;


public class Cell {
	private double startLat;
	private double endLat;
	private double startLng;
	private double endLng;
	private long startTime;
	private long endTime;
	private List<STPoint> points;
	private boolean clustered;
	
	public Cell(double startLat, double endLat, double startLng, double endLng, long startTime, long endTime) {
		this.startLat = startLat;
		this.endLat = endLat;
		this.startLng = startLng;
		this.endLng = endLng;
		this.startTime = startTime;
		this.endTime = endTime;
		points = new ArrayList<STPoint>();
		clustered = false;
	}
	
	public void addPoint(STPoint point) {
		points.add(point);
	}
	
	public List<STPoint> getPoints() {
		return points;
	}
	
	public int getPointsCount() {
		return points.size();
	}

	public double getStartLat() {
		return startLat;
	}

	public double getEndLat() {
		return endLat;
	}

	public double getStartLng() {
		return startLng;
	}

	public double getEndLng() {
		return endLng;
	}
	
	public long getStartTime() {
		return startTime;
	}
	
	public long getEndTime() {
		return endTime;
	}
	
	public boolean isClustered() {
		return clustered;
	}
	
	public void setClustered(boolean clustered) {
		this.clustered = clustered;
	}
	
	public String toString() {
		String description = "lat: " + getStartLat() + " .. " + getEndLat();
		description += " | lng: " + getStartLng() + " .. " + getEndLng();
		description += " | time: " + getStartTime() + " .. " + getEndTime();
		return description;
	}
}
