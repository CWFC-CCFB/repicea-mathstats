/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2022 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.model.glm.measerr;

import java.security.InvalidParameterException;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.data.DataSet;
import repicea.stats.data.GenericStatisticalDataStructure;
import repicea.stats.data.StatisticalDataStructure;
import repicea.stats.integral.TrapezoidalRule;
import repicea.stats.model.CompositeLogLikelihoodWithExplanatoryVariables;
import repicea.stats.model.IndividualLogLikelihood;
import repicea.stats.model.glm.LinkFunction;
import repicea.stats.model.glm.LinkFunction.Type;

/**
 * A class implementing a Berkson error model under the assumption of uniform distribution. <br>
 * <br>
 * This class relies on the trapezoidal rule for numerical integration.
 * 
 * @author Mathieu Fortin - July 2022
 */
abstract class GLMUniformBerksonMeasErrorDefinition extends AbstractGLMMeasErrorDefinition {

 
	private static class GLMWithUniformMeasErrorDataStructure extends GenericStatisticalDataStructure {

		private GLMWithUniformMeasErrorDataStructure(DataSet dataSet) {
			super(dataSet);
		}

		Matrix setMeasErrorDefinition(GLMUniformBerksonMeasErrorDefinition errDef) {
			Matrix errBounds = new Matrix(getNumberOfObservations(), 2);
			if (errDef.lowerBoundVar != null) {
				if (!dataSet.getFieldNames().contains(errDef.lowerBoundVar)) {
					throw new InvalidParameterException("The lowerBoundVar argument in the measurement error definition is not in the dataset!");
				} else {
					int indexField = dataSet.getIndexOfThisField(errDef.lowerBoundVar);
					for (int i = 0; i < dataSet.getNumberOfObservations(); i++) {
						errBounds.setValueAt(i, 0, (Double) dataSet.getObservations().get(i).getValueAt(indexField));
					}
				}
			}

			String upperBoundFieldName = errDef.upperBoundVar != null ? errDef.upperBoundVar : errDef.effectWithMeasError;
			if (!dataSet.getFieldNames().contains(upperBoundFieldName)) {
				throw new InvalidParameterException("The upperBound argument in the measurement error definition is not in the dataset!");
			} else {
				int indexField = dataSet.getIndexOfThisField(upperBoundFieldName);
				for (int i = 0; i < dataSet.getNumberOfObservations(); i++) {
					double maxValue = (Double) dataSet.getObservations().get(i).getValueAt(indexField);
					if (maxValue < errBounds.getValueAt(i, 0))
						throw new InvalidParameterException("The upper bound argument seems to contain values lower than the lower bound!");
					errBounds.setValueAt(i, 1, maxValue);
				}
			}
			return errBounds;
		}
	}
	
	
	@SuppressWarnings("serial")
	private static class GLMWithUniformMeasErrorCompositeLogLikelihood extends CompositeLogLikelihoodWithExplanatoryVariables {
		
		private final GLMUniformBerksonMeasErrorDefinition measErr;
		private final LinkFunctionWithMeasError lf;
		
		private GLMWithUniformMeasErrorCompositeLogLikelihood(IndividualLogLikelihood innerLogLikelihoodFunction, 
				Matrix xValues, 
				Matrix yValues, 
				GLMUniformBerksonMeasErrorDefinition measErr,
				LinkFunctionWithMeasError lf) {
			super(innerLogLikelihoodFunction, xValues, yValues);
			this.measErr = measErr;
			this.lf = lf;
		}
		
		@Override
		protected void setValuesInLikelihoodFunction(int index) {
			super.setValuesInLikelihoodFunction(index);
			lf.setBounds(measErr.bounds.getValueAt(index, 0), measErr.bounds.getValueAt(index, 1));
		}
	}

	@SuppressWarnings("serial")
	private static class LinkFunctionWithMeasError extends LinkFunction {

		private final LinkFunction linkFunctionErrorFreeObs;
		private final GradientHessianProvider gradientProvider;
		private final GradientHessianProvider hessianProvider;
		private double lowerBound;
		private double upperBound;
		
		private final TrapezoidalRule adaptedTr;
		private final GLMUniformBerksonMeasErrorDefinition measErr;
		
		private LinkFunctionWithMeasError(Type linkFunctionType, double resolution, GLMUniformBerksonMeasErrorDefinition measErr) {
			super(linkFunctionType);
			linkFunctionErrorFreeObs = new LinkFunction(linkFunctionType, getOriginalFunction());
			gradientProvider = new GradientHessianProvider(linkFunctionErrorFreeObs, true);
			hessianProvider = new GradientHessianProvider(linkFunctionErrorFreeObs, false);
			adaptedTr = new TrapezoidalRule(resolution);
			this.measErr = measErr;
		}
		
		void setBounds(double lowerBound, double upperBound) {
			this.lowerBound = lowerBound;
			this.upperBound = upperBound;
		}
		
		private boolean isObservationErrorFree() {
			return getOriginalFunction().getVariableValue(measErr.indexEffectWithMeasError) <= measErr.valueForNotConsideringMeasurementError;
		}
		
		@Override
		public Double getValue() {
			if (isObservationErrorFree()) {
				return linkFunctionErrorFreeObs.getValue();
			} else {
				double invdensity = upperBound - lowerBound;
				adaptedTr.setLowerBound(lowerBound);
				adaptedTr.setUpperBound(upperBound);
				double marginalProbability = adaptedTr.getIntegralApproximation(linkFunctionErrorFreeObs, measErr.indexEffectWithMeasError, false) / invdensity; // false: it is a variable.
				return marginalProbability;
			}
		}

		@Override
		public Matrix getGradient() {
			if (isObservationErrorFree()) {
				return linkFunctionErrorFreeObs.getGradient();
			} else {
				double invdensity = upperBound - lowerBound;
				adaptedTr.setLowerBound(lowerBound);
				adaptedTr.setUpperBound(upperBound);
				Matrix marginalGradient = adaptedTr.getIntegralApproximationForMatrixFunction(gradientProvider, measErr.indexEffectWithMeasError, false).scalarMultiply(1d/invdensity); // false: it is a variable.
				return marginalGradient;
			}
		}

		@Override
		public SymmetricMatrix getHessian() {
			if (isObservationErrorFree()) {
				return linkFunctionErrorFreeObs.getHessian();
			} else {
				double invdensity = upperBound - lowerBound;
				adaptedTr.setLowerBound(lowerBound);
				adaptedTr.setUpperBound(upperBound);
				Matrix marginalHessian = adaptedTr.getIntegralApproximationForMatrixFunction(hessianProvider, measErr.indexEffectWithMeasError, false).scalarMultiply(1d/invdensity); // false: it is a variable.
				return SymmetricMatrix.convertToSymmetricIfPossible(marginalHessian);
			}
		}
	}

	private final double valueForNotConsideringMeasurementError;
	private final String lowerBoundVar;
	private final String upperBoundVar;
	private final double resolution;
	
	private Matrix bounds;
	
	/**
	 * Constructor.
	 * @param effectWithMeasError the field name of the variable with measurement error
	 * @param valueForNotConsideringMeasurementError a value for which the measurement error is assumed to be negligible
	 * @param lowerBoundVar a field name representing the minimum value for the measurement error
	 * @param upperBoundVar a field name representing the maximum value for the measurement error
	 * @param resolution the resolution of the trapezoidal numerical integration 
	 */
	public GLMUniformBerksonMeasErrorDefinition(String effectWithMeasError, 
			double valueForNotConsideringMeasurementError,
			String lowerBoundVar,
			String upperBoundVar,
			double resolution) {
		super(MeasurementErrorModel.Berkson, effectWithMeasError);
		if (resolution <= 0) {
			throw new InvalidParameterException("The resolution argument must be strictly positive!");
		} else {
			this.resolution = resolution;
		}
		this.valueForNotConsideringMeasurementError = valueForNotConsideringMeasurementError;
		this.lowerBoundVar = lowerBoundVar;
		this.upperBoundVar = upperBoundVar;
	}

	@Override
	public void validate(GLMWithMeasurementError glm) {
		super.validate(glm);
		bounds = ((GLMWithUniformMeasErrorDataStructure) glm.getDataStructure()).setMeasErrorDefinition(this);
	}

	@Override
	public CompositeLogLikelihoodWithExplanatoryVariables createCompositeLikelihoodFromModel(GLMWithMeasurementError glm) {
		return new GLMWithUniformMeasErrorCompositeLogLikelihood(glm.getIndividualLogLikelihood(), 
				glm.getMatrixX(), 
				glm.getVectorY(),
				this, 
				(LinkFunctionWithMeasError) glm.getLinkFunction());
	}

	@Override
	public StatisticalDataStructure createDataStructureFromDataSet(DataSet dataSet) {
		return new GLMWithUniformMeasErrorDataStructure(dataSet);
	}

	@Override
	public IndividualLogLikelihood createIndividualLogLikelihoodFromModel(GLMWithMeasurementError glm) {
		return null;	// the default super method will be used instead
	}

	@Override
	public LinkFunction createLinkFunction(Type linkFunctionType, GLMWithMeasurementError glm) {
		return new LinkFunctionWithMeasError(linkFunctionType, resolution, this);
	}

	
}
