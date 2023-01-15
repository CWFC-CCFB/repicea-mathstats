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

import repicea.stats.integral.CompositeSimpsonRule;
import repicea.stats.integral.AbstractNumericalIntegrationMethod;
import repicea.stats.integral.TrapezoidalRule;

/**
 * The StemTaperSegment class represents a segment in the tree. It
 * provides two heights: one for the bottom cross section and a second
 * one for the top cross section.
 * @author Mathieu Fortin - July 2012
 */
public class StemTaperSegment {

	public static final double VERY_SMALL = 1E-3;

	private AbstractNumericalIntegrationMethod nim;
	
	/**
	 * Constructor 1.
	 * @param bottomCrossSection a double that represents the height of the lower cross section (m)
	 * @param topCrossSection a double that represents the height of the upper cross section (m)
	 * @param nim a NumericalIntegrationMethod instance
	 */
	public StemTaperSegment(double bottomCrossSection, double topCrossSection, AbstractNumericalIntegrationMethod nim) {
		if (bottomCrossSection > 0 && bottomCrossSection < topCrossSection) {
			nim.setLowerBound(bottomCrossSection);
			nim.setUpperBound(topCrossSection);
			this.nim = nim;
		} else {
			throw new InvalidParameterException("The bottom cross section is higher than the top section or its height is negative!");
		}
	}

	/**
	 * Constructor 2. 
	 * @param nim a NumericalIntegrationMethod instance that already contains the lower bound and upper bounds
	 */
	public StemTaperSegment(AbstractNumericalIntegrationMethod nim) {
		if (nim.getLowerBound() > 0 && nim.getLowerBound() < nim.getUpperBound()) {
			this.nim = nim;
		} else {
			throw new InvalidParameterException("The bottom cross section is higher than the top section or its height is negative!");
		}
	}
	
	
	/**
	 * Constructor 2.
	 * @param bottonStemTaperCrossSection a StemTaperCrossSection instance that represents the lower cross section
	 * @param topStemTaperCrossSection a StemTaperCrossSection instance that represents the upper cross section
	 * @param nim a NumericalIntegrationMethod instance
	 */
	protected StemTaperSegment(StemTaperCrossSection bottonStemTaperCrossSection, 
			StemTaperCrossSection topStemTaperCrossSection,
			AbstractNumericalIntegrationMethod nim) {
		this(bottonStemTaperCrossSection.getSectionHeight(), topStemTaperCrossSection.getSectionHeight(), nim);
	}

	/**
	 * This method returns the height (m) of the lower cross section.
	 * @return a double
	 */
	protected double getBottomHeight() {return nim.getLowerBound();}
	
	/**
	 * This method returns the height (m) of the upper cross section.
	 * @return a double
	 */
	protected double getTopHeight() {return nim.getUpperBound();}

	/**
	 * This method returns a List of Double that corresponds to the height for which
	 * the stem taper must be predicted. It depends on the integration method.
	 * @return a List of Double instances
	 */
	protected List<Double> getHeightsToEvaluate() {
		return nim.getXValues();
	}
	
	protected List<Double> getWeights() {
		return nim.getWeights();
	}
	
	protected List<Double> getRescalingFactors() {
		return nim.getRescalingFactors();
	}
	
	/**
	 * This method creates a StemTaperSegment with either a TrapezoidalRule or a CompositeSimpsonRule instance
	 * for numerical integration method.
	 * @param heights the list of heights in the segment (a List of Double instances)
	 * @param simpson true to use the CompositeSimpsonRule
	 * @param upToIndex the index in the heights to be considered
	 * @return a StemTaperSegment
	 */
	static StemTaperSegment createSegment(List<Double> heights, boolean simpson, int upToIndex) {
		if (simpson && upToIndex % 2 != 1) {
			upToIndex--;
		}
		List<Double> dumpList = new ArrayList<Double>();
		for (int i = 0; i < upToIndex; i++) {
			dumpList.add(heights.get(i));
		}
		StemTaperSegment segment;
		if (simpson) {
			segment = new StemTaperSegment(new CompositeSimpsonRule(dumpList));
		} else {
			segment = new StemTaperSegment(new TrapezoidalRule(dumpList));
		}
		dumpList.remove(dumpList.size() - 1);			// do not remove last element as it will serve for the next segment
		heights.removeAll(dumpList);
		return segment;
	}
	

	
//	public static void main(String[] args) {
//		List<Double> heights = new ArrayList<Double>();
//		heights.add(0.3);
//		heights.add(0.7);
//		heights.add(1.7);
//		heights.add(2.7);
//		heights.add(3.7);
//		heights.add(4.2);
//		StemTaperSegment.createStemTaperSegmentList(heights);
//	}
}
