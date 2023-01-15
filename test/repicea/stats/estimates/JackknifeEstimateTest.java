package repicea.stats.estimates;

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.util.ObjectUtility;

public class JackknifeEstimateTest {

	@SuppressWarnings("unchecked")
	@Test
	public void testDelete1Jackknife() {
		List<Double> observations = new ArrayList<Double>();
		for (double i = 1; i <= 10; i++) {
			observations.add(i);
		}

		Matrix mat = new Matrix(observations);
		double mean = mat.getSumOfElements() / mat.m_iRows;
		Matrix diff = mat.scalarAdd(-mean);
		double sumSquaredDiff = diff.transpose().multiply(diff).getSumOfElements();
		double estimatedVariance = sumSquaredDiff / ((mat.m_iRows - 1) * mat.m_iRows);
		
		Matrix tmpMat;
		Matrix realization;
		JackknifeEstimate estimate = new JackknifeEstimate(10, 1);
		for (int i = 1; i <= 10; i++) {
			List<Double> clonedList = ObjectUtility.copyList(observations);
			clonedList.remove(i - 1);
			tmpMat = new Matrix(clonedList);
			realization = new Matrix(1,1);
			realization.setValueAt(0, 0, tmpMat.getSumOfElements() / tmpMat.m_iRows);
			estimate.addRealization(realization);
		}
		
		double variance = estimate.getVariance().getValueAt(0, 0);
		Assert.assertEquals("Testing variance", estimatedVariance, variance, 1E-8);
	}

	
	
	@SuppressWarnings("unchecked")
	@Test
	public void testDelete2Jackknife() {
		List<Double> observations = new ArrayList<Double>();
		for (double i = 1; i <= 10; i++) {
			observations.add(i);
		}

		Matrix mat = new Matrix(observations);
		double mean = mat.getSumOfElements() / mat.m_iRows;
		Matrix diff = mat.scalarAdd(-mean);
		double sumSquaredDiff = diff.transpose().multiply(diff).getSumOfElements();
		double estimatedVariance = sumSquaredDiff / ((mat.m_iRows - 1) * mat.m_iRows);
		
		Matrix tmpMat;
		Matrix realization;
		JackknifeEstimate estimate = new JackknifeEstimate(10, 2);
		for (Double i = 1d; i <= (10 - 1); i++) {
			for (Double j = i+1; j <= 10; j++) {
				List<Double> clonedList = ObjectUtility.copyList(observations);
				clonedList.remove(i);
				clonedList.remove(j);
				tmpMat = new Matrix(clonedList);
				realization = new Matrix(1,1);
				realization.setValueAt(0, 0, tmpMat.getSumOfElements() / tmpMat.m_iRows);
				estimate.addRealization(realization);
			}
		}
		
		double variance = estimate.getVariance().getValueAt(0, 0);
		Assert.assertEquals("Testing variance", estimatedVariance, variance, 1E-8);
	}

	
	
	@SuppressWarnings("unchecked")
	@Test
	public void testDelete3Jackknife() {
		List<Double> observations = new ArrayList<Double>();
		for (double i = 1; i <= 10; i++) {
			observations.add(i);
		}

		Matrix mat = new Matrix(observations);
		double mean = mat.getSumOfElements() / mat.m_iRows;
		Matrix diff = mat.scalarAdd(-mean);
		double sumSquaredDiff = diff.transpose().multiply(diff).getSumOfElements();
		double estimatedVariance = sumSquaredDiff / ((mat.m_iRows - 1) * mat.m_iRows);
		
		Matrix tmpMat;
		Matrix realization;
		JackknifeEstimate estimate = new JackknifeEstimate(10, 3);
		for (Double i = 1d; i <= (10 - 2); i++) {
			for (Double j = i+1; j <= (10 - 1); j++) {
				for (Double k = j+1; k <= 10; k++) {
					List<Double> clonedList = ObjectUtility.copyList(observations);
					clonedList.remove(i);
					clonedList.remove(j);
					clonedList.remove(k);
					tmpMat = new Matrix(clonedList);
					realization = new Matrix(1,1);
					realization.setValueAt(0, 0, tmpMat.getSumOfElements() / tmpMat.m_iRows);
					estimate.addRealization(realization);
				}
			}
		}
		
		double variance = estimate.getVariance().getValueAt(0, 0);
		Assert.assertEquals("Testing variance", estimatedVariance, variance, 1E-8);
	}

}
