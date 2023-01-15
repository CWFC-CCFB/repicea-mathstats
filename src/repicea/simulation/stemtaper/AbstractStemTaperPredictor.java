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

import java.util.List;

import repicea.simulation.REpiceaPredictor;

/**
 * This interface is the basis for any stem taper model.
 * @author Mathieu Fortin - April 2014
 */
@SuppressWarnings("serial")
public abstract class AbstractStemTaperPredictor extends REpiceaPredictor {

	protected AbstractStemTaperPredictor(boolean isParametersVariabilityEnabled, boolean isRandomEffectsVariabilityEnabled,	boolean isResidualVariabilityEnabled) {
		super(isParametersVariabilityEnabled, isRandomEffectsVariabilityEnabled, isResidualVariabilityEnabled);
	}


	public abstract interface BasicStemTaperTree {}
	
	/**
	 * This method computes the stem taper.
	 * @param tree a BasicStemTaperTrees
	 * @param heightMeasures a List of Double that represent the height (m)
	 * @param additionalParameters a series of Object that are additional parameters
	 * @return an Estimate instance
	 */
	public abstract AbstractStemTaperEstimate getPredictedTaperForTheseHeights(BasicStemTaperTree tree, List<Double> heightMeasures, Object... additionalParameters);

	
	/**
	 * This method computes the stem taper.
	 * @param stemTaperSegments a List of StemTaperSegment instances
	 * @return a StemTaperEstimate instance with the cross section diameter
	 */
	public AbstractStemTaperEstimate getPredictedTaperForTheseSegments(BasicStemTaperTree tree, StemTaperSegmentList stemTaperSegments) {		
		List<Double> currentHeightsToEvaluate = stemTaperSegments.getHeightsWithoutReplicates();	
		return getPredictedTaperForTheseHeights(tree, currentHeightsToEvaluate);		
	}

	
}
