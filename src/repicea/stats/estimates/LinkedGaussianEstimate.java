package repicea.stats.estimates;

import java.util.ArrayList;
import java.util.List;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;

@SuppressWarnings("serial")
public class LinkedGaussianEstimate extends GaussianEstimate {

	private final List<Integer> indexFirstEstimate;
	private final List<Integer> indexSecondEstimate;
	
	public LinkedGaussianEstimate(GaussianEstimate firstEstimate, GaussianEstimate secondEstimate, double correlation) {
		super();
		indexFirstEstimate = new ArrayList<Integer>();
		indexSecondEstimate = new ArrayList<Integer>();
		int nbParms1 = firstEstimate.getMean().m_iRows;
		int nbParms2 = secondEstimate.getMean().m_iRows;
		for (int i = 0; i < (nbParms1 + nbParms2); i++) {
			if (i < nbParms1) {
				indexFirstEstimate.add(i);
			} else {
				indexSecondEstimate.add(i);
			}
		}
		
		Matrix mean = firstEstimate.getMean().matrixStack(secondEstimate.getMean(), true);
		SymmetricMatrix variance = new SymmetricMatrix(mean.m_iRows);
		variance.setSubMatrix(firstEstimate.getVariance(), 0, 0);
		variance.setSubMatrix(secondEstimate.getVariance(), nbParms1, nbParms1);
		for (int i = nbParms1; i < (nbParms1 + nbParms2); i++) {
			for (int j = 0; j < nbParms1; j++) {
				double var1 = variance.getValueAt(i, i);
				double var2 = variance.getValueAt(j, j);
				variance.setValueAt(i, j, Math.sqrt(var1 * var2) * correlation);
//				variance.setValueAt(j, i, variance.getValueAt(i, j));
			}
		}
		this.setMean(mean);
		this.setVariance(variance);
	}
	
	
//	public static void main(String[] args) {
//		Matrix mean1 = new Matrix(2,1);
//		mean1.m_afData[0][0] = 1d;
//		mean1.m_afData[1][0] = 2d;
//
//		Matrix variance1 = new Matrix(2,2);
//		variance1.m_afData[0][0] = 1d;
//		variance1.m_afData[1][0] = .5;
//		variance1.m_afData[0][1] = .5;
//		variance1.m_afData[1][1] = .75;
//
//		Matrix mean2 = new Matrix(3,1);
//		mean2.m_afData[0][0] = 3d;
//		mean2.m_afData[1][0] = 4d;
//		mean2.m_afData[2][0] = 5d;
//
//		Matrix variance2 = new Matrix(3,3);
//		variance2.m_afData[0][0] = .5;
//		variance2.m_afData[1][0] = .25;
//		variance2.m_afData[2][0] = .1;
//		variance2.m_afData[0][1] = .25;
//		variance2.m_afData[1][1] = 2d;
//		variance2.m_afData[2][1] = .15;
//		variance2.m_afData[0][2] = .1;
//		variance2.m_afData[1][2] = .15;
//		variance2.m_afData[2][2] = 1d;
//
//		LinkedGaussianEstimate estimate = new LinkedGaussianEstimate(new GaussianEstimate(mean1, variance1), new GaussianEstimate(mean2, variance2), .4);
//		
//	}
}
