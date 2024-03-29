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
 * @author Mathieu Fortin - March 2012
 * @param <M> an AbstractMatrix-derived class that stands for the mean
 * @param <V> an AbstractMatrix-derived class that stands for the variance
 * @param <D> a Distribution derived instance which represents the assumed distribution for the estimate
 */
@SuppressWarnings("rawtypes")
public abstract class Estimate<M extends AbstractMatrix, V extends AbstractMatrix, D extends Distribution<M,V>> extends RandomVariable<M,V,D> {
	
	static {
		SerializerChangeMonitor.registerEnumNameChange("repicea.stats.estimates.Estimate$EstimatorType", "MonteCarlo", "Resampling");
	}
	
	
	private static final long serialVersionUID = 20120825L;
	
	protected EstimatorType estimatorType;

	protected final List<String> rowIndex;
	
	/**
	 * The type of estimator.
	 * @author Mathieu Fortin - March 2012
s	 */
	public static enum EstimatorType {
		Resampling, 
		LeastSquares, 
		LikelihoodBased, 
		MomentBased, 
		Unknown}

	protected Estimate(D distribution) {
		super(distribution);
		rowIndex = new ArrayList<String>();
	}
	
	
	/**
	 * This method returns the type of the estimator.
	 * @return an EstimatorType instance
	 */
	public EstimatorType getEstimatorType() {return estimatorType;}

	/**
	 * This method makes it possible to set an optional row index. This is useful when the 
	 * response is a vector. 
	 * @param newRowIndex a List of String instance. A null value reset the 
	 */
	public void setRowIndex(List<String> newRowIndex) {
		this.rowIndex.clear();
		if (newRowIndex != null && !newRowIndex.isEmpty()) {
			if (newRowIndex.size() != getMean().m_iRows) {
				throw new InvalidParameterException("The size of the list is incompatible with tne dimension of the estimate!");
			}
			this.rowIndex.addAll(newRowIndex);
		}
	}

	
	/**
	 * This method returns a copy of the row index. 
	 * @return a List of String instance or null if the row index has not been set.
	 */
	public List<String> getRowIndex() {
		List<String> rowIndexCopy = new ArrayList<String>();
		rowIndexCopy.addAll(rowIndex);
		return rowIndexCopy;
	}
	
	/**
	 * This method returns a random deviate from this estimate. This method
	 * is useful for Monte Carlo simulations.
	 * @return a deviate from the underlying distribution as a Matrix instance
	 */
	public M getRandomDeviate() {
		return getDistribution().getRandomRealization();
	}

	
	/**
	 * This method returns an difference estimate.
	 * @param estimate2 an Estimate to be subtracted from this estimate.
	 * @return an Estimate
	 */
	@SuppressWarnings("unchecked")
	public Estimate<M, V, ?> getDifferenceEstimate(Estimate<M, V, ?> estimate2) {
		M diff = (M) getMean().subtract(estimate2.getMean());
		if (diff instanceof ComplexMatrix) {
			throw new UnsupportedOperationException("The method getDifferenceEstimate has not been implemented for ComplexMatrix yet!");
		}
		Matrix variance = (Matrix) getVariance().add(estimate2.getVariance());
		if (variance instanceof SymmetricMatrix) {
			return (Estimate<M, V, ?>) new GaussianEstimate((Matrix) diff, (SymmetricMatrix) variance);
		} else {
			throw new UnsupportedOperationException("The variance object is not a SymmetricMatrix instance!");
		}
	}
	
	/**
	 * This method returns a sum of two estimates.
	 * @param estimate2 an Estimate to be added to this estimate.
	 * @return an Estimate
	 */
	@SuppressWarnings("unchecked")
	public Estimate<M, V, ?> getSumEstimate(Estimate<M, V, ?> estimate2) {
		M diff = (M) getMean().add(estimate2.getMean());
		if (diff instanceof ComplexMatrix) {
			throw new UnsupportedOperationException("The method getSumEstimate has not been implemented for ComplexMatrix yet!");
		}
		Matrix variance = (Matrix) getVariance().add(estimate2.getVariance());
		if (variance instanceof SymmetricMatrix) {
			return (Estimate<M, V, ?>) new GaussianEstimate((Matrix) diff, (SymmetricMatrix) variance);
		} else {
			throw new UnsupportedOperationException("The variance object is not a SymmetricMatrix instance!");
		}
			
	}

	
	/**
	 * This method returns an difference estimate.
	 * @param scalar a double to be multiplied by this estimate
	 * @return an Estimate
	 */
	@SuppressWarnings("unchecked")
	public Estimate<M, V, ?> getProductEstimate(double scalar) {
		M diff = (M) getMean().scalarMultiply(scalar);
		if (diff instanceof ComplexMatrix) {
			throw new UnsupportedOperationException("The method getProductEstimate has not been implemented for ComplexMatrix yet!");
		}
		SymmetricMatrix variance = (SymmetricMatrix) getVariance().scalarMultiply(scalar * scalar);
		return (Estimate<M, V, ?>) new GaussianEstimate((Matrix) diff, variance);
	}

	/**
	 * This method returns the probability of getting a lower valueand upper bound of a confidence intervals at probability
	 * level 1 - alpha
	 * @param oneMinusAlpha is 1 minus the probability of Type I error
	 * @return a ConfidenceInterval instance 
	 */
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
	

	/**
	 * Returns an estimate of the product of two parametric univariate estimate. The variance
	 * estimator is based on Goodman's estimator.
	 * @param estimate an Estimate instance
	 * @return a SimpleEstimate instance
	 */
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
	
	/**
	 * Collapse the estimate following a map that contains the indices for each group. <br>
	 * <br>
	 * The collapsing ensures the consistency, that is all the row indices must be found in the
	 * list instances contained in the map argument. If there is a mismatch, the method will throw an 
	 * exception. IMPORTANT: the new indices, that is the keys of the map argument are sorted in the
	 * new Estimate instance.
	 * @param desiredIndicesForCollapsing a LinkedHashMap with the keys being the new indices and 
	 * the values being lists of indices to be collapsed.
	 * @return an Estimate instance
	 */
	public Estimate<Matrix, SymmetricMatrix, ?> collapseEstimate(LinkedHashMap<String, List<String>> desiredIndicesForCollapsing) {
		return collapseMeanAndVariance(desiredIndicesForCollapsing);
	}
	
	protected final Estimate<Matrix, SymmetricMatrix, ?> collapseMeanAndVariance(LinkedHashMap<String, List<String>> desiredIndicesForCollapsing) {
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
		
		Estimate<Matrix, SymmetricMatrix, ?> outputEstimate = new SimpleEstimate(newMean, SymmetricMatrix.convertToSymmetricIfPossible(newVariance));
		
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
