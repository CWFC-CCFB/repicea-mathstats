package repicea.stats.estimates;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.stats.Distribution.Type;
import repicea.stats.StatisticalUtility;

public class MultivariateAndCollapseTest {

	@Test
	public void completeCollapsingWithMonteCarloEstimateTest() {
		int nElements = 10;
		MonteCarloEstimate est = new MonteCarloEstimate();
		for (int i = 0; i < 1000; i++) {
			est.addRealization(StatisticalUtility.drawRandomVector(nElements, Type.GAUSSIAN));
		}
		
		Matrix basicMean = est.getMean();
		Assert.assertTrue("Testing we have 10 rows", basicMean.m_iRows == nElements);
		Matrix basicVariance = est.getVariance();
		Assert.assertTrue("Testing we have 10 rows", basicVariance.m_iRows == nElements);
		Assert.assertTrue("Testing we have 10 cols", basicVariance.m_iCols == nElements);
		
		List<String> rowIndex = new ArrayList<String>();
		for (int i = 0; i < nElements; i++) {
			rowIndex.add("" + i);
		}
		est.setRowIndex(rowIndex);
		LinkedHashMap<String, List<String>> collapseIndices = new LinkedHashMap<String, List<String>>();
		collapseIndices.put("all", rowIndex);
		Estimate<?> collapsedEstimate = est.collapseEstimate(collapseIndices);
		
		Matrix collapsedMean = collapsedEstimate.getMean();
		Assert.assertTrue("Testing we have 1 row", collapsedMean.m_iRows == 1);
		Assert.assertEquals("Testing consistency", basicMean.getSumOfElements(), collapsedMean.getValueAt(0, 0), 1E-8);
		Matrix collapsedVariance = collapsedEstimate.getVariance();
		Assert.assertTrue("Testing we have 1 row", collapsedVariance.m_iRows == 1);
		Assert.assertTrue("Testing we have 1 col", collapsedVariance.m_iCols == 1);
		Assert.assertEquals("Testing consistency", basicVariance.getSumOfElements(), collapsedVariance.getValueAt(0, 0), 1E-8);
	}

	@Test
	public void partialCollapsingWithMonteCarloEstimateTest() {
		int nElements = 10;
		MonteCarloEstimate est = new MonteCarloEstimate();
		for (int i = 0; i < 1000; i++) {
			est.addRealization(StatisticalUtility.drawRandomVector(nElements, Type.GAUSSIAN));
		}
		
		Matrix basicMean = est.getMean();
		Assert.assertTrue("Testing we have 10 rows", basicMean.m_iRows == nElements);
		Matrix basicVariance = est.getVariance();
		Assert.assertTrue("Testing we have 10 rows", basicVariance.m_iRows == nElements);
		Assert.assertTrue("Testing we have 10 cols", basicVariance.m_iCols == nElements);
		
		List<String> rowIndex = new ArrayList<String>();
		for (int i = 0; i < nElements; i++) {
			rowIndex.add("" + i);
		}
		est.setRowIndex(rowIndex);
		LinkedHashMap<String, List<String>> collapseIndices = new LinkedHashMap<String, List<String>>();
		collapseIndices.put("group1", new ArrayList<String>());
		collapseIndices.put("group2", new ArrayList<String>());
		for (int i = 0; i < rowIndex.size(); i++) {
			if (i < 3) {
				collapseIndices.get("group1").add(rowIndex.get(i));
			} else {
				collapseIndices.get("group2").add(rowIndex.get(i));
			}
		}
			
		Estimate<?> collapsedEstimate = est.collapseEstimate(collapseIndices);
		
		Matrix collapsedMean = collapsedEstimate.getMean();
		Assert.assertTrue("Testing we have 2 rows", collapsedMean.m_iRows == 2);
		Assert.assertEquals("Testing consistency", basicMean.getSumOfElements(), collapsedMean.getSumOfElements(), 1E-8);
		Matrix collapsedVariance = collapsedEstimate.getVariance();
		Assert.assertTrue("Testing we have 2 rows", collapsedVariance.m_iRows == 2);
		Assert.assertTrue("Testing we have 2 cols", collapsedVariance.m_iCols == 2);
		Assert.assertEquals("Testing consistency", basicVariance.getSumOfElements(), collapsedVariance.getSumOfElements(), 1E-8);
	}

}
