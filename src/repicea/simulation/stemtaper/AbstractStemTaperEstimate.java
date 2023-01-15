/*
 * This file is part of the repicea-foresttools library.
 *
 * Copyright (C) 2009-2014 Mathieu Fortin for Rouge-Epicea
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed with the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * Please see the license at http://www.gnu.org/copyleft/lesser.html.
 */
package repicea.simulation.stemtaper;

import java.security.InvalidParameterException;
import java.util.List;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.CentralMomentsSettable;
import repicea.stats.Distribution;
import repicea.stats.distributions.GaussianDistribution;
import repicea.stats.estimates.Estimate;
import repicea.stats.estimates.GaussianEstimate;
import repicea.stats.integral.TrapezoidalRule;

/**
 * The StemTaperEstimate class extends the Estimate class. It handles MonteCarlo or maximum likelihood based estimates.
 * likelihood based. The method getMean() and getVariance() are overriden in order to select which super methods apply.
 * The StemTaperEstimate handles the volume integration.
 * @author Mathieu Fortin - March 2012
 * @author Mathieu Fortin - December 2016 (Refactoring)
 */
@SuppressWarnings({ "rawtypes", "serial" })
public abstract class AbstractStemTaperEstimate extends Estimate<Distribution> implements CentralMomentsSettable  {
		
	private List<Double> heights;
	private Estimate<?> volumeEstimate;


	/**
	 * Constructor.
	 * @param computedHeights a List instance containing the heights (m) of the cross sections
	 */
	public AbstractStemTaperEstimate(List<Double> computedHeights) {
		super(new GaussianDistribution(null, null));
		estimatorType = EstimatorType.LikelihoodBased;
		heights = computedHeights;
	}

	
	/**
	 * This method post processes the predictions. If a particular height appears twice in the originalHeightsToEvaluate
	 * member, the matrix row (and column) corresponding to this heights is doubled. 
	 * @param originalPredictions the original predictions (Matrix instance)
	 * @return the post processing predictions (also Matrix instance)
	 */
	private Matrix reshapeMatrixAccordingToSegments(StemTaperSegmentList segments, Matrix originalPredictions) {
		List<Double> originalHeightsToEvaluate = segments.getHeights();
		
		if (originalHeightsToEvaluate.size() == heights.size() && heights.containsAll(originalHeightsToEvaluate)) {
			return originalPredictions;
		} else {
			Matrix outputMatrix;
			boolean isColumnVector = false;
			if (originalPredictions.isColumnVector()) {
				isColumnVector = true;
				outputMatrix = new Matrix(originalHeightsToEvaluate.size(), 1);
			} else if (originalPredictions.isSquare()) {
				outputMatrix = new Matrix(originalHeightsToEvaluate.size(), originalHeightsToEvaluate.size());
			} else {
				throw new InvalidParameterException("Matrix originalPredictions is not square nor a row vector!");
			}
			int i_index, j_index;
			for (int i = 0; i < outputMatrix.m_iRows; i++) {
				i_index = heights.indexOf(originalHeightsToEvaluate.get(i));
				if (isColumnVector) {
					j_index = 0;
					outputMatrix.setValueAt(i, 0, originalPredictions.getValueAt(i_index, j_index));
				} else {		// is square then
					for (int j = i; j < outputMatrix.m_iCols; j++) {
						j_index = heights.indexOf(originalHeightsToEvaluate.get(j));
						outputMatrix.setValueAt(i, j, originalPredictions.getValueAt(i_index, j_index));
						outputMatrix.setValueAt(j, i, originalPredictions.getValueAt(i_index, j_index));
					}
				}
			}
			return outputMatrix;
		}
	}
	

	/**
	 * This method returns a volume estimate from the integrated taper all along the heights contained
	 * in the heights member.
	 * @return an Estimate instance 
	 */
	public Estimate<?> getVolumeEstimate() {
		return getVolumeEstimate(null);
	}
	
	
	/**
	 * This method returns a volume estimate from the integrated taper all along the selected segments.
	 * @param segments a StemTaperSegmentList that represents the selected taper (can be null, in such case the integration is carried out
	 * all along the heights, by default)
	 * @return an Estimate instance 
	 */
	public Estimate<?> getVolumeEstimate(StemTaperSegmentList segments) {
		if (segments == null) {	// means it has to be integrated all along the sections
			if (volumeEstimate == null) {
				segments = getDefaultSegments();
				volumeEstimate = getVolumeEstimateForTheseSegments(segments);
			} 
			return volumeEstimate;
		} else {
			return getVolumeEstimateForTheseSegments(segments);
		}
	}

	
	private StemTaperSegmentList getDefaultSegments() {
		StemTaperSegmentList segments = new StemTaperSegmentList();
		segments.add(new StemTaperSegment(new TrapezoidalRule(heights)));
		return segments;
	}


	/**
	 * This method returns a volume estimate from the integrated taper all along the selected segments.
	 * @param segments a StemTaperSegmentList instance
	 * @return an Estimate instance 
	 */
	private Estimate<?> getVolumeEstimateForTheseSegments(StemTaperSegmentList segments) {
		if (!heights.containsAll(segments.getHeights())) {
			throw new InvalidParameterException("There is a mismatch between the requested sections and the heights that have been computed!");
		}

		Matrix weights = new Matrix(segments.getWeightsAcrossSegments());
		Matrix rescalingFactors = new Matrix(segments.getRescalingFactorsAcrossSegments());
		Matrix volumeFactor = rescalingFactors.elementWiseMultiply(weights).scalarMultiply(getScalingFactor());
		Matrix varianceFactor = volumeFactor.multiply(volumeFactor.transpose());

		Estimate result = new GaussianEstimate();
		((GaussianEstimate) result).setMean(null);
		((GaussianEstimate) result).setVariance(null);

		SymmetricMatrix variance;

		Matrix taper = getSquaredDiameters(reshapeMatrixAccordingToSegments(segments, getMean()));

		Matrix volumeEstim = taper.elementWiseMultiply(volumeFactor);
		((GaussianEstimate) result).setMean(volumeEstim);
		if (getVariance() != null) {
			variance = getVarianceOfSquaredDiameter(
					SymmetricMatrix.convertToSymmetricIfPossible(
							reshapeMatrixAccordingToSegments(segments, getVariance())));
			Matrix scaledVariance = variance.elementWiseMultiply(varianceFactor);
			if (scaledVariance instanceof SymmetricMatrix) {
				((GaussianEstimate) result).setVariance((SymmetricMatrix) scaledVariance);
			} else {
				throw new UnsupportedOperationException("The scaledVariance object is not a SymmetricMatrix instance!");
			}
		}
		return result;
	}


	/**
	 * This method returns the heights of the cross sections for which the square diameters were predicted.
	 * @return a List of Double instances
	 */
	public List<Double> getCrossSectionHeights() {return heights;}
	
	protected abstract Matrix getSquaredDiameters(Matrix predictedDiameters);
	
	protected abstract SymmetricMatrix getVarianceOfSquaredDiameter(SymmetricMatrix variancePredictedDiameters);
	
	protected abstract double getScalingFactor();
	
	@Override
	public void setMean(Matrix mean) {
		((GaussianDistribution) getDistribution()).setMean(mean);
	}

	@Override
	public void setVariance(SymmetricMatrix variance) {
		((GaussianDistribution) getDistribution()).setVariance(variance);
	}

}
