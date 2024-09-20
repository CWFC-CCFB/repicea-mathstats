/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin - Canadian Wood Fibre Centre, CFS
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This library is distributed with the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * Please see the license at http://www.gnu.org/copyleft/lesser.html.
 */
package repicea.stats.model.lm;

import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.stats.data.DataSet;
import repicea.stats.model.lm.LogBackTransformation.Estimator;
import repicea.util.ObjectUtility;

/**
 * Used for manuscript Fortin 2023. Log transformation and Regression.
 * @author Mathieu Fortin - August 2023
 *
 */
public class CaseStudyTruncatedGaussian {

	private static Object[] convertMatrixToObjectArray(Matrix m) {
		Object[] predFieldValues = new Object[m.m_iRows];
		for (int i = 0; i < m.m_iRows; i++)
			predFieldValues[i] = m.getValueAt(i, 0);
		return predFieldValues;
	}
	
	@Test
	public void predictionTest() throws Exception {
		String inputFilename = ObjectUtility.getPackagePath(CaseStudyTruncatedGaussian.class) + "datasetSingleObs.csv";
		DataSet ds = new DataSet(inputFilename, true);
		String formula = "yTrans ~ lnDt_corr + dbhCm + BAL";
		LinearModel lm = new LinearModel(ds, formula);
		lm.doEstimation();
		ds.addField("meanPredOLS", convertMatrixToObjectArray(lm.getPredicted()));
		Matrix predOriginal = LogBackTransformation.getMeanPredictedValuesOnOriginalScale(lm, 1, Estimator.Naive);
		ds.addField("meanPredOLSOrigScale", convertMatrixToObjectArray(predOriginal));
		System.out.println(lm.getSummary());
		Matrix p = lm.getParameters();
		Matrix parmsTrad = new Matrix(p.m_iRows + 1, 1);
		parmsTrad.setSubMatrix(p, 0, 0);
		parmsTrad.setValueAt(parmsTrad.m_iRows - 1, 0, lm.getResidualVariance());
		LinearModelWithTruncatedGaussianErrorTerm truncatedModel = new LinearModelWithTruncatedGaussianErrorTerm(ds, formula, parmsTrad, 0);
		truncatedModel.doEstimation();
		ds.addField("meanPredMML", convertMatrixToObjectArray(truncatedModel.getPredicted()));
		predOriginal = LogBackTransformation.getMeanPredictedValuesOnOriginalScale(truncatedModel, 1, Estimator.Naive);
		ds.addField("meanPredMMLOrigScale", convertMatrixToObjectArray(predOriginal)); // offset = 1 and no variance
		System.out.println(truncatedModel.getSummary());
		
		String outputDatasetFilename = ObjectUtility.getPackagePath(CaseStudyTruncatedGaussian.class) + "datasetSingleObsOutput.csv";
		DataSet ref_ds = new DataSet(outputDatasetFilename, true);
		
		List<Object> actual = ds.getFieldValues(ds.getIndexOfThisField("meanPredMMLOrigScale"));
		List<Object> expected = ref_ds.getFieldValues(ds.getIndexOfThisField("meanPredMMLOrigScale"));
		Assert.assertEquals("Testing nb of observations", expected.size(), actual.size());
		for (int i = 0; i < actual.size(); i++) {
			Assert.assertEquals("Testing values", (Double) expected.get(i), (Double) actual.get(i), 0.1);
		}
//		String outputDatasetFilename = ObjectUtility.getPackagePath(CaseStudyTruncatedGaussian.class) + "datasetSingleObsOutput.csv";
//		outputDatasetFilename = outputDatasetFilename.replace("/bin", "");
//		ds.save(outputDatasetFilename);
	}
	
	
}
