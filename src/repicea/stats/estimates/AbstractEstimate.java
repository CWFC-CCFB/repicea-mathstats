/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2012 Mathieu Fortin for Rouge-Epicea
 * Copyright (C) 2025 His Majesty the King in right of Canada
 * Author: Mathieu Fortin, Canadian Forest Service
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
package repicea.stats.estimates;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;

import repicea.math.AbstractMatrix;
import repicea.math.ComplexMatrix;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.serial.SerializerChangeMonitor;
import repicea.stats.Distribution;
import repicea.stats.RandomVariable;

/**
 * The Estimate class is the basic class for all estimates.
 * @author Mathieu Fortin - March 2012, January 2025
 * @param <M> an AbstractMatrix-derived class that stands for the mean
 * @param <V> an AbstractMatrix-derived class that stands for the variance
 * @param <D> a Distribution derived instance which represents the assumed distribution for the estimate
 */
@SuppressWarnings({ "rawtypes", "serial" })
public abstract class AbstractEstimate<M extends AbstractMatrix, V extends AbstractMatrix, D extends Distribution<M,V>> extends RandomVariable<M,V,D>
									implements Estimate<M,V,D> {
	
	static {
		SerializerChangeMonitor.registerEnumNameChange("repicea.stats.estimates.Estimate$EstimatorType", "MonteCarlo", "Resampling");
	}
	
//	private static final long serialVersionUID = 20120825L;
	
	protected EstimatorType estimatorType;

	protected final List<String> rowIndex;
	

	protected AbstractEstimate(D distribution) {
		super(distribution);
		rowIndex = new ArrayList<String>();
	}
	
	
	@Override
	public EstimatorType getEstimatorType() {return estimatorType;}

	@Override
	public void setRowIndex(List<String> newRowIndex) {
		this.rowIndex.clear();
		if (newRowIndex != null && !newRowIndex.isEmpty()) {
			if (newRowIndex.size() != getMean().m_iRows) {
				throw new InvalidParameterException("The size of the list is incompatible with tne dimension of the estimate!");
			}
			this.rowIndex.addAll(newRowIndex);
		}
	}

	
	@Override
	public List<String> getRowIndex() {
		List<String> rowIndexCopy = new ArrayList<String>();
		rowIndexCopy.addAll(rowIndex);
		return rowIndexCopy;
	}
	
	@Override
	public M getRandomDeviate() {
		return getDistribution().getRandomRealization();
	}

	
	@SuppressWarnings("unchecked")
	@Override
	public Estimate<M, V, ?> getDifferenceEstimate(Estimate<M, V, ?> estimate2) {
		M diff = (M) getMean().subtract(estimate2.getMean());
		if (diff instanceof ComplexMatrix) {
			throw new UnsupportedOperationException("The method getDifferenceEstimate has not been implemented for ComplexMatrix yet!");
		}
		Matrix variance = (Matrix) getVariance().add(estimate2.getVariance());
		if (variance instanceof SymmetricMatrix) {
			return (AbstractEstimate<M, V, ?>) new GaussianEstimate((Matrix) diff, (SymmetricMatrix) variance);
		} else {
			throw new UnsupportedOperationException("The variance object is not a SymmetricMatrix instance!");
		}
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public Estimate<M, V, ?> getSumEstimate(Estimate<M, V, ?> estimate2) {
		M diff = (M) getMean().add(estimate2.getMean());
		if (diff instanceof ComplexMatrix) {
			throw new UnsupportedOperationException("The method getSumEstimate has not been implemented for ComplexMatrix yet!");
		}
		Matrix variance = (Matrix) getVariance().add(estimate2.getVariance());
		if (variance instanceof SymmetricMatrix) {
			return (AbstractEstimate<M, V, ?>) new GaussianEstimate((Matrix) diff, (SymmetricMatrix) variance);
		} else {
			throw new UnsupportedOperationException("The variance object is not a SymmetricMatrix instance!");
		}
			
	}

	@SuppressWarnings("unchecked")
	@Override
	public Estimate<M, V, ?> getProductEstimate(double scalar) {
		M diff = (M) getMean().scalarMultiply(scalar);
		if (diff instanceof ComplexMatrix) {
			throw new UnsupportedOperationException("The method getProductEstimate has not been implemented for ComplexMatrix yet!");
		}
		SymmetricMatrix variance = (SymmetricMatrix) getVariance().scalarMultiply(scalar * scalar);
		return (AbstractEstimate<M, V, ?>) new GaussianEstimate((Matrix) diff, variance);
	}

	@Override
	public abstract ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha);
	
	/**
	 * This method checks if the two point estimates are compatible. The basic
	 * check consists of comparing the classes. Then, the matrix data is checked
	 * for consistency with previous data.
	 * @param estimate an Estimate instance
	 * @return a boolean
	 */
	protected boolean isMergeableEstimate(Estimate<?,?,?> estimate) {
		return false;
	}
	

	@Override
	public SimpleEstimate getProductEstimate(Estimate<M, V, ?> estimate) {
		if (estimate.getDistribution().isUnivariate() && getDistribution().isUnivariate()) {
			M alphaMu = getMean();
			M betaMu = estimate.getMean();
			if (alphaMu instanceof ComplexMatrix || betaMu instanceof ComplexMatrix) {
				throw new UnsupportedOperationException("The method getProductEstimate has not been implemented for ComplexMatrix yet!");
			}
			Matrix alphaMean = (Matrix) alphaMu;
			Matrix betaMean = (Matrix) betaMu;
			Matrix alphaVariance = (Matrix) getVariance();
			Matrix betaVariance = (Matrix) estimate.getVariance();
			Matrix newMean = alphaMean.multiply(betaMean);
			Matrix newVariance = alphaMean.elementWisePower(2d).multiply(betaVariance).
					add(betaMean.elementWisePower(2d).multiply(alphaVariance)).
					subtract(alphaVariance.multiply(betaVariance));
			return new SimpleEstimate(newMean,  SymmetricMatrix.convertToSymmetricIfPossible(newVariance));
		}
		throw new InvalidParameterException("The getProductEstimate is only implemented for parametric univariate distribution ");
	}
	
	/**
	 * A static method to compute the product of many estimates.
	 * @param estimates
	 * @return a SimpleEstimate instance
	 */
	public static SimpleEstimate getProductOfManyEstimates(List<Estimate<Matrix, SymmetricMatrix, ?>> estimates) {
		Estimate<Matrix, SymmetricMatrix, ?> currentEstimate = null;
		for (int i = 1; i < estimates.size(); i++) {
			if (i == 1) {
				currentEstimate = estimates.get(i-1);
			} 
			currentEstimate = currentEstimate.getProductEstimate(estimates.get(i));
		}
		return (SimpleEstimate) currentEstimate;
	}
	
	@Override
	public AbstractEstimate<Matrix, SymmetricMatrix, ?> collapseEstimate(LinkedHashMap<String, List<String>> desiredIndicesForCollapsing) {
		return collapseMeanAndVariance(desiredIndicesForCollapsing);
	}
	
	protected final AbstractEstimate<Matrix, SymmetricMatrix, ?> collapseMeanAndVariance(LinkedHashMap<String, List<String>> desiredIndicesForCollapsing) {
		M m = getMean();
		if (!(m instanceof Matrix)) {
			throw new UnsupportedOperationException("The collapseEstimate method is meant for Estimate based on Matrix instance!");
		}
		Matrix mean = (Matrix) m;
		if (rowIndex.isEmpty()) {
			throw new InvalidParameterException("The row indices have not been set yet!");
		}
		if (rowIndex.size() != mean.m_iRows) {
			throw new InvalidParameterException("The size of the list is incompatible with tne dimension of the estimate!");
		}
		List<String> copyOfIndex = new ArrayList<String>();
		copyOfIndex.addAll(getRowIndex());
		Collections.sort(copyOfIndex);
		List<String> completeList = new ArrayList<String>();
		for (List<String> l : desiredIndicesForCollapsing.values()) {
			completeList.addAll(l);
		}
		Collections.sort(completeList);
		if (!completeList.equals(copyOfIndex)) {
			throw new InvalidParameterException("Some indices are missing in the desiredIndicesForCollapsing or cannot be found in the row indices!");
		} 

//		P oldMean = getMean();
		Matrix newMean = collapseRowVector(mean, desiredIndicesForCollapsing);
		
		Matrix oldVariance = (Matrix) getVariance();
		Matrix newVariance = collapseSquareMatrix(oldVariance, desiredIndicesForCollapsing);
		
		AbstractEstimate<Matrix, SymmetricMatrix, ?> outputEstimate = new SimpleEstimate(newMean, SymmetricMatrix.convertToSymmetricIfPossible(newVariance));
		
		List<String> newIndexRow = new ArrayList<String>(desiredIndicesForCollapsing.keySet());
		Collections.sort(newIndexRow);
		outputEstimate.setRowIndex(newIndexRow);
		
		return outputEstimate;
	}

	
	protected final Matrix collapseRowVector(Matrix originalMatrix, LinkedHashMap<String, List<String>> desiredIndicesForCollapsing) {
		List<String> newIndexRow = new ArrayList<String>(desiredIndicesForCollapsing.keySet());
		Collections.sort(newIndexRow);
		Matrix collapsedMatrix = new Matrix(desiredIndicesForCollapsing.size(), 1);
		for (int i = 0; i < collapsedMatrix.m_iRows; i++) {
			List<String> requestedIndices = desiredIndicesForCollapsing.get(newIndexRow.get(i));
			collapsedMatrix.setValueAt(i, 0, originalMatrix.getSubMatrix(convertIndexIntoInteger(requestedIndices), null).getSumOfElements());
		}
		return collapsedMatrix;
	}

	protected final Matrix collapseSquareMatrix(Matrix originalMatrix, LinkedHashMap<String, List<String>> desiredIndicesForCollapsing) {
		if (originalMatrix == null) { // means the variance cannot be calculated as in the MonteCarloEstimate class if the nb realizations is smaller than 2.
			return null;
		} else {
			List<String> newIndexRow = new ArrayList<String>(desiredIndicesForCollapsing.keySet());
			Collections.sort(newIndexRow);
			Matrix collapsedMatrix = new Matrix(desiredIndicesForCollapsing.size(), desiredIndicesForCollapsing.size());
			for (int i = 0; i < collapsedMatrix.m_iRows; i++) {
				List<String> requestedIndices_i = desiredIndicesForCollapsing.get(newIndexRow.get(i));
				for (int j = 0; j < collapsedMatrix.m_iRows; j++) {
					List<String> requestedIndices_j = desiredIndicesForCollapsing.get(newIndexRow.get(j));
					Matrix tmpMatrix = originalMatrix.getSubMatrix(convertIndexIntoInteger(requestedIndices_i), 
							convertIndexIntoInteger(requestedIndices_j));
					collapsedMatrix.setValueAt(i, j, tmpMatrix.getSumOfElements());
				}
			}
			return collapsedMatrix;
		}
	}
	
	
	private List<Integer> convertIndexIntoInteger(List<String> selectedIndices) {
		List<Integer> outputList = new ArrayList<Integer>();
		for (String s : selectedIndices) {
			outputList.add(getRowIndex().indexOf(s));
		}
		return outputList;
	}
	
	@Override
	public final M getMean() {
		return super.getMean();
	}

	@Override
	public final V getVariance() {
		return super.getVariance();
	}
	
}
