import ij.*;

public class FilteringSession {

	/*******************************************************************************
	 *
	 * E D G E   D E T E C T O R   S E C T I O N
	 *
	 ******************************************************************************/

	/**
	 * Detects the vertical edges inside an ImageAccess object.
	 * This is the non-separable version of the edge detector.
	 * The kernel of the filter has the following form:
	 *
	 *     -------------------
	 *     | -1  |  0  |  1  |
	 *     -------------------
	 *     | -1  |  0  |  1  |
	 *     -------------------
	 *     | -1  |  0  |  1  |
	 *     -------------------
	 *
	 * Mirror border conditions are applied.
	 */
	static public ImageAccess detectEdgeVertical_NonSeparable(ImageAccess input) {
		int nx = input.getWidth();
		int ny = input.getHeight();
		double arr[][] = new double[3][3];
		double pixel;
		ImageAccess out = new ImageAccess(nx, ny);
		for (int x = 0; x < nx; x++) {
			for (int y = 0; y < ny; y++) {
				input.getNeighborhood(x, y, arr);
				pixel = arr[2][0]+arr[2][1]+arr[2][2]-arr[0][0]-arr[0][1]-arr[0][2];
				pixel = pixel / 6.0;
				out.putPixel(x, y, pixel);
			}
		}
		return out;
	}

	/**
	 * Detects the vertical edges inside an ImageAccess object.
	 * This is the separable version of the edge detector.
	 * The kernel of the filter applied to the rows has the following form:
	 *     -------------------
	 *     | -1  |  0  |  1  |
	 *     -------------------
	 *
	 * The kernel of the filter applied to the columns has the following 
	 * form:
	 *     -------
	 *     |  1  |
	 *     -------
	 *     |  1  |
	 *     -------
	 *     |  1  |
	 *     -------
	 *
	 * Mirror border conditions are applied.
	 */
	static public ImageAccess detectEdgeVertical_Separable(ImageAccess input) {
		int nx = input.getWidth();
		int ny = input.getHeight();
		ImageAccess out = new ImageAccess(nx, ny);
		double rowin[]  = new double[nx];
		double rowout[] = new double[nx];
		for (int y = 0; y < ny; y++) {
			input.getRow(y, rowin);
			doDifference3(rowin, rowout);
			out.putRow(y, rowout);
		}
		
		double colin[]  = new double[ny];
		double colout[] = new double[ny];
		for (int x = 0; x < nx; x++) {
			out.getColumn(x, colin);
			doAverage3(colin, colout);
			out.putColumn(x, colout);
		}
		return out;
	}

	static public ImageAccess detectEdgeHorizontal_NonSeparable(ImageAccess input) {
		int nx = input.getWidth();
		int ny = input.getHeight();
		double arr[][] = new double[3][3];
		double pixel;
		ImageAccess out = new ImageAccess(nx, ny);
		for (int x = 0; x < nx; x++) {
			for (int y = 0; y < ny; y++) {
				input.getNeighborhood(x, y, arr);
				pixel = -arr[0][0]-arr[1][0]-arr[2][0]+arr[0][2]+arr[1][2]+arr[2][2];
				pixel = pixel / 6.0;
				out.putPixel(x, y, pixel);
			}
		}
		return out;
	}
	

	static public ImageAccess detectEdgeHorizontal_Separable(ImageAccess input) {
		int nx = input.getWidth();
		int ny = input.getHeight();
		ImageAccess out = new ImageAccess(nx, ny);
		double rowin[]  = new double[nx];
		double rowout[] = new double[nx];
		for (int y = 0; y < ny; y++) {
			input.getRow(y, rowin);
			doAverage3(rowin, rowout);
			out.putRow(y, rowout);
		}
		
		double colin[]  = new double[ny];
		double colout[] = new double[ny];
		for (int x = 0; x < nx; x++) {
			out.getColumn(x, colin);
			doDifference3(colin, colout);
			out.putColumn(x, colout);
		}
		return out;

	}

	/**
	 * Implements an one-dimensional average filter of length 3.
	 * The filtered value of a pixel is the averaged value of
	 * its local neighborhood of length 3.
	 * Mirror border conditions are applied.
	 */
	static private void doAverage3(double vin[], double vout[]) {
		int n = vin.length;
		vout[0] = (vin[0] + 2.0 * vin[1]) / 3.0;
		for (int k = 1; k < n-1; k++) {
			vout[k] = (vin[k-1] + vin[k] + vin[k+1]) / 3.0;
		}
		vout[n-1] = (vin[n-1] + 2.0 * vin[n-2]) / 3.0;
	}

	/**
	 * Implements an one-dimensional centered difference filter of 
	 * length 3. The filtered value of a pixel is the difference of 
	 * its two neighborhing values.
	 * Mirror border conditions are applied.
	 */
	static private void doDifference3(double vin[], double vout[]) {
		int n = vin.length;
		vout[0] = 0.0;
		for (int k = 1; k < n-1; k++) {
			vout[k] = (vin[k+1] - vin[k-1]) / 2.0;
		}
		vout[n-1] = 0.0;
	}
	
	static private void doAverage5(double vin[], double vout[]) {
		int n = vin.length;
		vout[0] = (vin[0] + 2.0 * vin[1] + 3.0 * vin[2]) / 6.0;
		vout[1] = (vin[1] + 2.0 * vin[2] + 3.0 * vin[3]) / 6.0;
		for(int k = 2; k < n-2; k++) {
			vout[k] = (vin[k-2]+ vin[k-1] + vin[k] + vin[k+1] + vin[k+2]) / 5.0;
		}	
		vout[n-2] = (vin[n-2] + 2.0 * vin[n-3] + 3.0 * vin[n-4]) / 6.0;
		vout[n-1] = (vin[n-1] + 2.0 * vin[n-2] + 3.0 * vin[n-3]) / 6.0;
	}
	
	static private void doAverage5_recursive(double vin[], double vout[]) {
		int n = vin.length;		
		vout[0] = (vin[0] + 2.0 * vin[1] + 3.0 * vin[2]) / 6.0;
		vout[1] = (vin[1] + 2.0 * vin[2] + 3.0 * vin[3]) / 6.0;
		for(int k = 2; k < n-2; k++) {
			vout[k] = soma_recursive(vin, k+2, 4)/5.0; 
		}
		vout[n-1] = (vin[n-1] + 2.0 * vin[n-2] + 3.0 * vin[n-3]) / 6.0;
		vout[n-2] = (vin[n-2] + 2.0 * vin[n-3] + 3.0 * vin[n-4]) / 6.0;					
	}
	
	static private double soma_recursive(double vin[],int z, int y) {
		double soma = 0.0; 
		if (y == 0){
			return soma;	
		}
		else{
			soma = vin[z] + soma_recursive(vin, z-1, y-1);
		}
		return soma;
	}
	/*******************************************************************************
	 *
	 * M O V I N G   A V E R A G E   5 * 5   S E C T I O N
	 *
	 ******************************************************************************/

	static public ImageAccess doMovingAverage5_NonSeparable(ImageAccess input) {
		int nx = input.getWidth();
		int ny = input.getHeight();
		ImageAccess out = new ImageAccess(nx, ny);	
		double arr[][] = new double[5][5];
		double pixel = 0.0;
			for(int x = 0; x < nx ; x++)
			for(int y = 0; y < ny; y++){
				input.getNeighborhood(x, y, arr);
				for(int m = 0; m < 5; m++)
				for(int n = 0; n < 5; n++) {
					pixel = pixel + arr[m][n];
				}
				pixel = pixel / 25.0;
				out.putPixel(x, y, pixel);
			}
		return out;
	}

	static public ImageAccess doMovingAverage5_Separable(ImageAccess input) {
		int nx = input.getWidth();
		int ny = input.getHeight();
		ImageAccess out = new ImageAccess(nx, ny);
		double rowin[] = new double[nx];
		double rowout[] = new double[nx];
		for(int y = 0; y < ny; y++) {
			input.getRow(y, rowin);
			doAverage5(rowin, rowout);
			out.putRow(y, rowout);
		}
		double colin[] = new double[ny];
		double colout[] = new double[ny];
		for (int x = 0; x < nx; x++) {
			out.getColumn(x, colin);
			doAverage5(colin, colout);
			out.putColumn(x, colout);
		}
		return out;	
	}

	static public ImageAccess doMovingAverage5_Recursive(ImageAccess input) {
		int nx = input.getWidth();
		int ny = input.getHeight();
		ImageAccess out = new ImageAccess(nx, ny);
		double rowin[] = new double[nx];
		double rowout[] = new double[nx];
		for (int y = 0; y < ny; y++) {
			input.getRow(y, rowin);
			doAverage5_recursive(rowin, rowout);
			out.putRow(y, rowout);
		}
		double colin[] = new double[ny];
		double colout[] = new double[ny];
		for (int x = 0; x < nx; x++) {
			out.getColumn(x, colin);
			doAverage5_recursive(colin, colout);
			out.putColumn(x, colout);
		}
		return out;
	}

	/*******************************************************************************
	 *
	 * S O B E L
	 *
	 ******************************************************************************/

	static public ImageAccess doSobel(ImageAccess input) {
		int nx = input.getWidth();
		int ny = input.getHeight();
		double pixel1 = 0.0;
		double pixel2 = 0.0;
		double arr[][] = new double[3][3];
		ImageAccess out1 = new ImageAccess(nx, ny);
		ImageAccess out2 = new ImageAccess(nx, ny);
		for (int y = 0; y < ny; y++)
		for (int x = 0; x < nx; x++){
				input.getNeighborhood(x, y, arr);
				pixel1 = arr[2][0] + (2 * arr[2][1]) +arr[2][2] - arr[0][0] - (2 * arr[0][1]) - arr[0][2];
				pixel1 = pixel1/6.0;
				out1.putPixel(x, y, pixel1);
		}
		for (int x = 0; x < nx; x++)
		for (int y = 0; y < ny; y++) {
				input.getNeighborhood(x, y, arr);
				pixel2 = -arr[0][0] - (2 * arr[1][0]) - arr[2][0] + arr[0][2] + (2 * arr[1][2]) + arr[2][2];
				pixel2 = pixel2/6.0;
				out2.putPixel(x, y, pixel2);
		}	
		out1.pow(2.0);
		out2.pow(2.0);
		out2.add(out1,out2);
		out2.sqrt();
		return out2;
	}

	/*******************************************************************************
	 *
	 * M O V I N G   A V E R A G E   L * L   S E C T I O N
	 *
	 ******************************************************************************/

	static public ImageAccess doMovingAverageL_Recursive(ImageAccess input, int length) {
		IJ.showMessage("Question 5");
		return input.duplicate();
	}

}
