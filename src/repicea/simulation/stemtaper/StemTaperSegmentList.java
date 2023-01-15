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
import java.util.ArrayList;
import java.util.List;

public class StemTaperSegmentList extends ArrayList<StemTaperSegment> {

	private static final long serialVersionUID = 20120713L;

	/**
	 * This method returns a List of Double instance that corresponds to the heights (m) across all
	 * segments.
	 * @return a List of Double instances
	 */
	protected List<Double> getHeights() {
		List<Double> heights= new ArrayList<Double>();
		for (StemTaperSegment segments : this) {
			heights.addAll(segments.getHeightsToEvaluate());
		}

		return heights;
	}
	
	
	/**
	 * This method returns a List of Double instance that corresponds to the heights (m) across all
	 * segments without replicates.	This method is essential to avoid having the matrices singular.
	 * @return a List of Double instances
	 */
	public List<Double> getHeightsWithoutReplicates() {
		List<Double> heightsWithoutReplicates = new ArrayList<Double>();
		List<Double> originalHeights = getHeights();
		for (Double height : originalHeights) {
			if (!heightsWithoutReplicates.contains(height)) {
				heightsWithoutReplicates.add(height);
			}
		}
		return heightsWithoutReplicates;
	}

	
	/**
	 * This method returns a List of Double instances that correspond to the weights across the segments.
	 * @return a List of Double
	 */
	protected List<Double> getWeightsAcrossSegments() {
		List<Double> weights = new ArrayList<Double>();
		for (StemTaperSegment segment : this) {
			weights.addAll(segment.getWeights());
		}
		return weights;
	}

	
	/**
	 * This method returns a List of Double instances that correspond to the rescaling factors across the segments.
	 * @return a List of Double
	 */
	protected List<Double> getRescalingFactorsAcrossSegments() {
		List<Double> rescalingFactors = new ArrayList<Double>();
		for (StemTaperSegment segment : this) {
			rescalingFactors.addAll(segment.getRescalingFactors());
		}
		return rescalingFactors;
	}

	/**
	 * This static methods provides a list of StemTaperSegment instances embedded in a StemTaperSegmentList.
	 * The method optimizes the integration by using CompositeSimpsonRule numerical integration whenever possible. 
	 * Otherwise, it uses TrapezoidalRule numerical integration.
	 * @param heights a List of Double instances that represent the heights along the bole
	 * @param optimize if false it disables the CompositeSimpsonRule and only the TrapezoidalRule is enabled
	 * @return a StemTaperSegmentList instance
	 */
	public static StemTaperSegmentList createStemTaperSegmentList(List<Double> heights, boolean optimize) {
		if (heights == null || heights.size() <= 1) {
			throw new InvalidParameterException("The height list must contains at least two values!");
		}

		StemTaperSegmentList segments = new StemTaperSegmentList();
		List<Double> tempList = new ArrayList<Double>();
		boolean wereTwoLastEqual = false;
		boolean areTwoLastEqual = false;
		int numberOfHeightsToConsider;
		for (Double point : heights) {
			tempList.add(point);
			if (optimize) {
				if (tempList.size() > 2) {
					areTwoLastEqual = areDistancesEqual(tempList);
					if (tempList.size() == 3) {
						wereTwoLastEqual = areTwoLastEqual;
					}
					if (areTwoLastEqual && !wereTwoLastEqual) {			// means we have a trapezoidal rule here
						numberOfHeightsToConsider = tempList.size() - 2;
						segments.add(StemTaperSegment.createSegment(tempList, false, numberOfHeightsToConsider));
						areTwoLastEqual = true;
						wereTwoLastEqual = true;

					} else if (!areTwoLastEqual && wereTwoLastEqual) {
						numberOfHeightsToConsider = tempList.size() - 1;
						segments.add(StemTaperSegment.createSegment(tempList, true, numberOfHeightsToConsider));
						areTwoLastEqual = false;
						wereTwoLastEqual = false;
					}
				}
			}
		}
		while (tempList.size() > 1) {
			segments.add(StemTaperSegment.createSegment(tempList, wereTwoLastEqual, tempList.size()));
			wereTwoLastEqual = false;
		}
		return segments;
	}

	
	private static boolean areDistancesEqual(List<Double> points) {
		double previousDiff = Double.NaN;
		double diff;
		for (int i = points.size() - 2; i < points.size(); i++) {
			diff = points.get(i) - points.get(i - 1); 
			if (Double.isNaN(previousDiff)) {
				previousDiff = diff;
			} else {
				if (Math.abs(diff-previousDiff) > StemTaperSegment.VERY_SMALL) {
					return false;
				}
			}
		}
		return true;
	}

}
