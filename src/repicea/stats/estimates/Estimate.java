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

import java.util.LinkedHashMap;
import java.util.List;

import repicea.math.AbstractMatrix;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.MomentGettable;
import repicea.stats.Distribution;

/**
 * An interface that ensures that the Estimate instance can provide
 * the appropriate information.
 * @author Mathieu Fortin - March 2012, January 2025
 * @param <M> an AbstractMatrix-derived class that stands for the mean
 * @param <V> an AbstractMatrix-derived class that stands for the variance
 * @param <D> a Distribution derived instance which represents the assumed distribution for the estimate
 */
@SuppressWarnings("rawtypes")
public interface Estimate<M extends AbstractMatrix, V extends AbstractMatrix, D extends Distribution<M,V>> 
						extends MomentGettable<M,V>, DistributionProvider<D> {


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

	/**
	 * Provide the type of the estimator.
	 * @return an EstimatorType instance
	 */
	public EstimatorType getEstimatorType(); 

	/**
	 * Make it possible to set an optional row index.<p>
	 * This is useful when the response is a vector. 
	 * @param newRowIndex a List of String instance. A null value resets the row index.
	 */
	public void setRowIndex(List<String> newRowIndex);

	
	/**
	 * Return a copy of the row index. 
	 * @return a List of String instance or null if the row index has not been set.
	 */
	public List<String> getRowIndex();
	
	/**
	 * Provide a random deviate from this estimate.<p>
	 * This method is used in Monte Carlo simulations.
	 * @return a deviate from the underlying distribution as a Matrix instance
	 */
	public M getRandomDeviate();

	
	/**
	 * Return a difference estimate.
	 * @param estimate2 an Estimate to be subtracted from this estimate.
	 * @return an Estimate
	 */
	public Estimate<M, V, ?> getDifferenceEstimate(Estimate<M, V, ?> estimate2);
	
	/**
	 * Return a sum of two estimates.
	 * @param estimate2 an Estimate to be added to this estimate.
	 * @return an Estimate
	 */
	public Estimate<M, V, ?> getSumEstimate(Estimate<M, V, ?> estimate2);

	
	/**
	 * Return the product of this estimate by a scalar.
	 * @param scalar a double to be multiplied by this estimate
	 * @return an Estimate
	 */
	public Estimate<M, V, ?> getProductEstimate(double scalar);

	/**
	 * This method returns the probability of getting a lower valueand upper bound of a confidence intervals at probability
	 * level 1 - alpha
	 * @param oneMinusAlpha is 1 minus the probability of Type I error
	 * @return a ConfidenceInterval instance 
	 */
	public ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha);
	

	/**
	 * Returns an estimate of the product of two parametric univariate estimate. The variance
	 * estimator is based on Goodman's estimator.
	 * @param estimate an Estimate instance
	 * @return a SimpleEstimate instance
	 */
	public Estimate<M, V, ?> getProductEstimate(Estimate<M, V, ?> estimate);

	
	/**
	 * Collapse the estimate following a map that contains the indices for each group. <p>
	 * The collapsing ensures the consistency, that is all the row indices must be found in the
	 * list instances contained in the map argument. If there is a mismatch, the method will throw an 
	 * exception. IMPORTANT: the new indices, that is the keys of the map argument are sorted in the
	 * new Estimate instance.
	 * @param desiredIndicesForCollapsing a LinkedHashMap with the keys being the new indices and 
	 * the values being lists of indices to be collapsed.
	 * @return an Estimate instance
	 */
	public Estimate<Matrix, SymmetricMatrix, ?> collapseEstimate(LinkedHashMap<String, List<String>> desiredIndicesForCollapsing);
	
	

}
