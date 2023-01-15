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

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.BlockingDeque;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.logging.Level;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.data.DataSet;
import repicea.stats.estimates.Estimate;
import repicea.stats.estimates.GaussianEstimate;
import repicea.stats.estimates.MonteCarloEstimate;
import repicea.stats.estimators.AbstractEstimator;
import repicea.stats.estimators.AbstractEstimator.EstimatorCompatibleModel;
import repicea.stats.model.glm.measerr.SIMEXModel.InternalGLM;
import repicea.util.REpiceaLogManager;

class SIMEXEstimator extends AbstractEstimator<EstimatorCompatibleModel> {

	public static String LOGGER_NAME = "SIMEXEstimator";
	
	class InternalWorker extends Thread {
		
		final InternalGLM model;
		final Map<Double, MonteCarloEstimate> estimateMap;
		final Map<Double, MonteCarloEstimate> varianceMap;
		EstimatorException estimatorException;
		
		InternalWorker(int id, InternalGLM model, Map<Double, MonteCarloEstimate> estimateMap, Map<Double, MonteCarloEstimate> varianceMap) {
			super("InternalWorker no " + id);
			this.model = model;
			this.estimateMap = estimateMap;
			this.varianceMap = varianceMap;
			start();
		}

		@Override
		public void run() {
			try {
				while (!this.isInterrupted()) {
					Object o = queue.takeFirst();
					if (o.equals(finalToken)) {
						break;
					} else {
						double factor = (Double) o;
						model.getCompleteLogLikelihood().generateMeasurementError(factor);
						if (refParms != null) {
							model.setParameters(refParms);
						} 
						model.doEstimation();
						if (model.getEstimator().isConvergenceAchieved()) {
							addRealizationToEstimate(estimateMap, varianceMap, factor, model);
							if (factor == 0d) {
								refParms = model.getEstimator().getParameterEstimates().getMean();
							}
						} else {
							queue.addFirst(factor);		// else we put the factor back into the map to make sure that this realization 
														// is going to converge at some point MF20221004
						}
						
					}
				}
			} catch (Exception e) {
				estimatorException = new EstimatorException(e.getMessage());
				SIMEXEstimator.this.interruptTasks();
			}
		}
	}

	private final Object finalToken = new Object();

	private final BlockingDeque<Object> queue;
	private GaussianEstimate parameterEstimates;
	private boolean convergenceAchieved;
	private List<InternalWorker> threads;
	
	DataSet parmsObsDataSet;
	DataSet parmsPredDataSet;
	private Matrix refParms;
	
	protected SIMEXEstimator(EstimatorCompatibleModel model) {
		super(model);
		queue = new LinkedBlockingDeque<Object>();
	}

	void interruptTasks() {
		if (threads != null) {
			for (InternalWorker t : threads) {
				t.interrupt();
			}
		}
	}

	protected SIMEXModel getModel() {return (SIMEXModel) model;}
	
	synchronized void addRealizationToEstimate(Map<Double, MonteCarloEstimate> estimateMap, 
			Map<Double, MonteCarloEstimate> varianceMap,
			double factor, 
			InternalGLM model) {
		if (!estimateMap.containsKey(factor)) {
			estimateMap.put(factor, new MonteCarloEstimate());
		}
		estimateMap.get(factor).addRealization(model.getParameters());
		if (!varianceMap.containsKey(factor)) {
			varianceMap.put(factor, new MonteCarloEstimate());
		}
		varianceMap.get(factor).addRealization(model.getEstimator().getParameterEstimates().getVariance());
	}

	
	@Override
	public boolean doEstimation() throws EstimatorException {
		convergenceAchieved = false;
		parameterEstimates = null;
		Map<Double, MonteCarloEstimate> varianceMap = new HashMap<Double, MonteCarloEstimate>();
		Map<Double, MonteCarloEstimate> estimateMap = new HashMap<Double, MonteCarloEstimate>();
		threads = new ArrayList<InternalWorker>();
		for (int id = 1; id <= getModel().nbThreads; id++) {
			threads.add(new InternalWorker(id, getModel().originalGLM.clone(), estimateMap, varianceMap));
		}
		for (Double factor : getModel().factors) {
			for (int b = 0; b < getModel().nbBootstrapRealizations; b++) {
				queue.add(factor);	// just one realization for the 0 factor
				if (factor == 0) {
					break;
				}
			}
		}
		for (@SuppressWarnings("unused") InternalWorker t : threads) 
			queue.add(finalToken);
		for (InternalWorker t : threads)
			try {
				t.join();
			} catch (InterruptedException e) {}

		for (InternalWorker t : threads) {
			if (t.estimatorException != null) {
				REpiceaLogManager.logMessage(LOGGER_NAME, Level.SEVERE, LOGGER_NAME, t.estimatorException.getMessage());
				return false;
			}
		}

		Matrix epsilon = new Matrix(getModel().factors);
		// A quadratic extrapolation x = [1, epsilon, epsilon2]
		Matrix x = new Matrix(epsilon.m_iRows, 1, 1, 0).matrixStack(epsilon, false).matrixStack(epsilon.elementWisePower(2d), false);
		Matrix invXtX_Xt = x.transpose().multiply(x).getInverseMatrix().multiply(x.transpose());
		Matrix parameters = null;
		Matrix variances = null;
		for (Double d : getModel().factors) {
			Matrix theseParms = estimateMap.get(d).getMean().transpose();
			Matrix theseVCov;
			SymmetricMatrix varianceMapResult = SymmetricMatrix.convertToSymmetricIfPossible(varianceMap.get(d).getMean());
			if (d == 0) {
				theseVCov = varianceMapResult.symSquare().transpose(); 
			} else {
				Matrix s2_delta = estimateMap.get(d).getVariance();
				theseVCov = ((SymmetricMatrix) varianceMapResult.subtract(s2_delta)).symSquare().transpose(); 
			}
			if (parameters == null) {
				parameters = theseParms;
				variances = theseVCov;
			} else { 
				parameters = parameters.matrixStack(theseParms, true);
				variances = variances.matrixStack(theseVCov, true);
			}
		}
		
		this.getModel().getEffectList();
		
		//			System.out.println("Parameters = " + parameters.toString());
		Matrix simexParms = new Matrix(parameters.m_iCols, 1);
		parmsObsDataSet = new DataSet(Arrays.asList(new String[] {"parmID", "obs", "zeta"}));
		parmsPredDataSet = new DataSet(Arrays.asList(new String[] {"parmID", "pred", "zeta"}));
		int nbRows = (int) Math.ceil((getMaxFactor() + 1) / .01);
		Matrix epsilonPred = new Matrix(nbRows, 1, -1, .01);
		Matrix xPred = new Matrix(nbRows, 1, 1, 0).matrixStack(epsilonPred, false).matrixStack(epsilonPred.elementWisePower(2d), false);
		Object[] record;
		for (int j = 0; j < parameters.m_iCols; j++) {
			Matrix y = parameters.getSubMatrix(0, parameters.m_iRows - 1, j, j);
			for (int i = 0; i < parameters.m_iRows; i++) {
				record = new Object[3];
				int jj = model.isInterceptModel() ? j - 1 : j;
				String name;
				if (jj < model.getEffectList().size()) {
					name = jj == -1 ? "intercept" : model.getEffectList().get(jj);
				} else {
					name = model.getOtherParameterNames().get(jj - model.getEffectList().size());
				}
				record[0] = name;
				record[1] = y.getValueAt(i, 0);	// observed
				record[2] = epsilon.getValueAt(i, 0); // zeta value
				parmsObsDataSet.addObservation(record);
			}
			Matrix beta = invXtX_Xt.multiply(y);
			Matrix simexValue = xPred.multiply(beta);
			simexParms.setValueAt(j, 0, simexValue.getValueAt(0, 0));
			for (int i = 0; i < xPred.m_iRows; i++) {
				record = new Object[3];
				int jj = model.isInterceptModel() ? j - 1 : j;
				String name;
				if (jj < model.getEffectList().size()) {
					name = jj == -1 ? "intercept" : model.getEffectList().get(jj);
				} else {
					name = model.getOtherParameterNames().get(jj - model.getEffectList().size());
				}
				record[0] = name;
				record[1] = simexValue.getValueAt(i, 0);	// predicted
				record[2] = epsilonPred.getValueAt(i, 0); // zeta value
				parmsPredDataSet.addObservation(record);
			}
		}
		Matrix simexVCov = new Matrix(variances.m_iCols, 1);
		Matrix extrapolation = xPred.getSubMatrix(0, 0, 0, xPred.m_iCols - 1);
		for (int j = 0; j < variances.m_iCols; j++) {
			Matrix y = variances.getSubMatrix(0, variances.m_iRows - 1, j, j);
			Matrix beta = invXtX_Xt.multiply(y);
			Matrix simexValue = extrapolation.multiply(beta);
			simexVCov.setValueAt(j, 0, simexValue.getValueAt(0, 0));
		}
		SymmetricMatrix simexVCovSymm = simexVCov.squareSym();
		parameterEstimates = new GaussianEstimate(simexParms, simexVCovSymm);
		getModel().predGLM.getEstimator().parameterEstimates = parameterEstimates;
		convergenceAchieved = true;
		return true;
	}

	private double getMaxFactor() {
		double max = 0d;
		for (double d : this.getModel().factors) {
			if (d > max)
				max = d;
		}
		return max;
	}
	
	@Override
	public boolean isConvergenceAchieved() {return convergenceAchieved;}

	@Override
	public Estimate<?> getParameterEstimates() {return parameterEstimates;}

	
	@Override
	public String getReport() {
		if (!isConvergenceAchieved()) {
			return "Convergence could not be achieved!";
		} else {
			StringBuilder sb = new StringBuilder();
			DataSet convergenceDataset = getConvergenceStatusReport();
			DecimalFormat decFormat = new DecimalFormat();
			decFormat.setMaximumFractionDigits(4);
			decFormat.setMinimumFractionDigits(4);
			convergenceDataset.setFormatter(1, decFormat);
			sb.append(convergenceDataset.toString() + System.lineSeparator());
			DataSet parameterDataset = getParameterEstimatesReport();
			decFormat = new DecimalFormat();
			decFormat.setMaximumFractionDigits(6);
			decFormat.setMinimumFractionDigits(6);
			parameterDataset.setFormatter(1, decFormat);
			parameterDataset.setFormatter(2, decFormat);
			parameterDataset.setFormatter(4, decFormat);
			decFormat = new DecimalFormat();
			decFormat.setMaximumFractionDigits(3);
			decFormat.setMinimumFractionDigits(3);
			parameterDataset.setFormatter(3, decFormat);
			sb.append(parameterDataset.toString() + System.lineSeparator());
			return sb.toString();
		}

	}

	
}

