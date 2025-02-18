/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2012 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.estimators;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.data.DataSet;
import repicea.stats.estimates.AbstractEstimate;
import repicea.stats.estimates.GaussianEstimate;
import repicea.stats.estimates.VarianceEstimate;
import repicea.stats.estimators.AbstractEstimator.EstimatorCompatibleModel;
import repicea.stats.estimators.OLSEstimator.OLSCompatibleModel;
import repicea.stats.model.lm.LinearModel;

/**
 * The OLSOptimizer implements the Ordinary Least Squares estimator.
 * @author Mathieu Fortin - November 2012
 */
public final class OLSEstimator extends AbstractEstimator<OLSCompatibleModel> {

	public interface OLSCompatibleModel extends EstimatorCompatibleModel {
	
		public Matrix getMatrixX();
		
		public Matrix getVectorY();
		
		public Matrix getResiduals();
		
		public void setParameters(Matrix beta);
	}
	
	
	private VarianceEstimate residualVariance;
	private boolean hasConverged;
	private AbstractEstimate<Matrix, SymmetricMatrix, ?> betaVector;
	private SymmetricMatrix inverseProduct;
	
	/**
	 * Constructor.
	 * @param model an OLSCompatibleModel instance
	 */
	public OLSEstimator(OLSCompatibleModel model) {
		super(model);
	}
	
	@Override
	public boolean doEstimation() throws EstimatorException {
		hasConverged = false;
		if (!(model instanceof LinearModel)) {
			throw new EstimatorException("The OLS optimizer is designed to work with instances of LinearModel only!");
		}
//		dataStruct = model.getDataStructure();
		Matrix matrixX = model.getMatrixX();
		Matrix matrixY = model.getVectorY();
		Matrix matrixXT = matrixX.transpose();
		betaVector = new GaussianEstimate();
		inverseProduct = SymmetricMatrix.convertToSymmetricIfPossible(matrixXT.multiply(matrixX).getInverseMatrix());
		((GaussianEstimate) betaVector).setMean(inverseProduct.multiply(matrixX.transpose()).multiply(matrixY));
		model.setParameters(betaVector.getMean());
		hasConverged = true;
		Matrix residual = model.getResiduals();
		int degreesOfFreedom = model.getNumberOfObservations() - betaVector.getMean().m_iRows;
		double resVar = residual.transpose().multiply(residual).scalarMultiply(1d / degreesOfFreedom).getValueAt(0, 0);
		residualVariance = new VarianceEstimate(degreesOfFreedom, resVar);
		((GaussianEstimate) betaVector).setVariance(inverseProduct.scalarMultiply(resVar));
		return true;
	}

	/**
	 * Provides the residual variance.
	 * @return a VarianceEstimate instance
	 */
	public VarianceEstimate getResidualVariance() {
		return residualVariance;
	}

	/**
	 * Provide the inverse of the X^T X.
	 * @return a SymmetricMatrix instance
	 */
	public SymmetricMatrix getInverseXtXMatrix() {
		return inverseProduct;
	}
	
	@Override
	public boolean isConvergenceAchieved() {return hasConverged;}

	@Override
	public AbstractEstimate<Matrix, SymmetricMatrix, ?> getParameterEstimates() {
		return betaVector;
	}


	@Override
	public DataSet getConvergenceStatusReport() {
		NumberFormat formatter = NumberFormat.getInstance();
		formatter.setMaximumFractionDigits(3);
		formatter.setMinimumFractionDigits(3);
		List<String> fieldNames = new ArrayList<String>();
		fieldNames.add("Element");
		fieldNames.add("Value");
		DataSet dataSet = new DataSet(fieldNames);
		Object[] record = new Object[2];
		record[0] = "Converged";
		record[1] = isConvergenceAchieved();
		dataSet.addObservation(record);
		return dataSet;
	}

	public String getReport() {
		StringBuilder sb = new StringBuilder();
		sb.append(super.getReport());
		sb.append("Residual variance: " + residualVariance.getMean().getValueAt(0, 0) + " with " + residualVariance.getDegreesOfFreedom() + " d.f.");
		return sb.toString();
	}
	
}
