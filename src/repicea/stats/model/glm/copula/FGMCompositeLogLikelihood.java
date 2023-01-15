package repicea.stats.model.glm.copula;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import repicea.math.AbstractMathematicalFunctionWrapper;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.math.optimizer.NewtonRaphsonOptimizer;
import repicea.math.optimizer.OptimizerListener;
import repicea.math.utility.MatrixUtility;
import repicea.stats.data.DataBlock;
import repicea.stats.data.HierarchicalStatisticalDataStructure;
import repicea.stats.model.CompositeLogLikelihoodWithExplanatoryVariables;
import repicea.stats.model.IndividualLogLikelihood;
import repicea.stats.model.Likelihood;

@SuppressWarnings("serial")
public class FGMCompositeLogLikelihood extends CompositeLogLikelihoodWithExplanatoryVariables implements OptimizerListener {

	private final static double VERY_SMALL = 1E-8;

	protected double llk;
	protected boolean llkUptoDate;
	
	protected Map<List<Integer>, Double> additionalLlkTerm;
	protected boolean additionalLlkTermUptoDate;
	
	protected Matrix gradientVector;
	protected boolean gradientVectorUptoDate;
	
	protected Map<List<Integer>, Matrix> additionalGradients;
	protected boolean additionalGradientTermUptoDate;
	
	protected SymmetricMatrix hessianMatrix;
	protected boolean hessianMatrixUptoDate;
	
	protected Map<List<Integer>, SymmetricMatrix> additionalHessians;
	protected boolean additionalHessianTermUptoDate;

	protected HierarchicalStatisticalDataStructure hierarchicalStructure;
	protected CopulaExpression copulaExpression;
	
	protected FGMCompositeLogLikelihood(IndividualLogLikelihood innerLogLikelihoodFunction, 
			Matrix xValues, 
			Matrix yValues) {
		super(innerLogLikelihoodFunction, xValues, yValues);
	}
	
	protected void initialize(HierarchicalStatisticalDataStructure hierarchicalStructure, CopulaExpression copulaExpression) {
		this.hierarchicalStructure = hierarchicalStructure;
		this.copulaExpression = copulaExpression;
	}
		
	@Override
	public void reset() {
		llkUptoDate = false;
		additionalLlkTermUptoDate = false;
		gradientVectorUptoDate = false;
		additionalGradientTermUptoDate = false;
		hessianMatrixUptoDate = false;
		additionalHessianTermUptoDate = false;
	}
	

	private Likelihood getIndividualLikelihood() {return (Likelihood) ((AbstractMathematicalFunctionWrapper) getOriginalFunction()).getOriginalFunction();}
//	private int getTotalNumberOfParameters() {return getIndividualLikelihood().getNumberOfParameters() + copulaExpression.getNumberOfParameters();}
			
			
	@Override
	public Matrix getGradient() {
		if (!gradientVectorUptoDate) {
			Matrix gradient = new Matrix(getNumberOfParameters(), 1);
			gradient.setSubMatrix(super.getGradient(),0,0);		// get the gradient under the assumption of independence
			
			for (Matrix additionalGradient : getAdditionalGradients().values()) {
//				MatrixUtility.add(gradient, additionalGradient);			// get the additional part of the gradient on both the beta vector and the copula parameters
				gradient = gradient.add(additionalGradient);			// get the additional part of the gradient on both the beta vector and the copula parameters
			}
			
			gradientVector = gradient;
			gradientVectorUptoDate = true;
		}
		return gradientVector;
	}

	@Override
	public SymmetricMatrix getHessian() {
		if (!hessianMatrixUptoDate) {
			SymmetricMatrix hessian = new SymmetricMatrix(getNumberOfParameters());
			SymmetricMatrix originalHessian = super.getHessian();
			for (int i = 0; i < originalHessian.m_iRows; i++) {
				for (int j = i; j < originalHessian.m_iRows; j++) {
					hessian.setValueAt(i, j, originalHessian.getValueAt(i, j));
				}
			}
			
//			hessian.setSubMatrix(super.getHessian(), 0, 0); 	// get the hessian under the assumption of independence
			
			for (SymmetricMatrix additionalHessian : getAdditionalHessians().values()) {
//				MatrixUtility.add(hessian, additionalHessian);
				hessian = (SymmetricMatrix) hessian.add(additionalHessian);
			}
			
			hessianMatrix = hessian;
			hessianMatrixUptoDate = true;
		}
		return hessianMatrix;
	}

	@Override
	public Double getValue() {
		if (!llkUptoDate) {
			double logLikelihood = super.getValue();
			for (Double additionalTerm : getAdditionalLikelihoodTerm().values()) {
				logLikelihood += Math.log(additionalTerm);
			}
			llk = logLikelihood;
			llkUptoDate = true;
		}
		return llk;
	}
	
	
	private Map<List<Integer>, Double> getAdditionalLikelihoodTerm() {
		Map<List<Integer>, Double> results = new HashMap<List<Integer>, Double>();
		
		if (!additionalLlkTermUptoDate) {
			
			int indexFirstObservation;
			double likelihoodFirst;
			double observedFirst;
			
			int indexSecondObservation;
			double likelihoodSecond;
			double observedSecond;
			
			double sumObserved;
			double multiplyingFactor;
			
			Map<String, DataBlock> map = hierarchicalStructure.getHierarchicalStructure();
			for (DataBlock db : map.values()) {
				List<Integer> index = db.getIndices();
				double additionalTerm = 1d;

				for (int i = 0; i < index.size() - 1; i++) {
					indexFirstObservation = index.get(i);
					setValuesInLikelihoodFunction(indexFirstObservation);
					likelihoodFirst = getIndividualLikelihood().getValue();
					observedFirst =  getIndividualLikelihood().getYVector().getValueAt(0, 0);
					for (int j = i + 1; j < index.size(); j++) {
						indexSecondObservation = index.get(j);
						setValuesInLikelihoodFunction(indexSecondObservation);
						likelihoodSecond =  getIndividualLikelihood().getValue();
						observedSecond =  getIndividualLikelihood().getYVector().getValueAt(0, 0);

						sumObserved = observedFirst + observedSecond;

						multiplyingFactor = 1d;
						if (Math.abs(sumObserved - 1) < VERY_SMALL) {
							multiplyingFactor = -1d;
						}

						copulaExpression.setX(indexFirstObservation, indexSecondObservation);
						double copulaValue = copulaExpression.getValue();

						additionalTerm += multiplyingFactor * copulaValue * (1 - likelihoodFirst) * (1 - likelihoodSecond);
					}
				}
				results.put(index, additionalTerm);
			}
			additionalLlkTerm = results;
			additionalLlkTermUptoDate = true;
		}
		return additionalLlkTerm;
	}
	
	
	private Map<List<Integer>, Matrix> getAdditionalGradients() {
		if (!additionalGradientTermUptoDate) {
			Map<List<Integer>, Matrix> additionalGradients = new HashMap<List<Integer>, Matrix>();
			
			int indexFirstObservation;
			double likelihoodFirst;
			Matrix du_dbetaFirst;
			double observedFirst;
			
			int indexSecondObservation;
			double likelihoodSecond;
			Matrix du_dbetaSecond;
			double observedSecond;
			
			double sumObserved;
			double multiplyingFactor;

			Matrix tmp;

			Map<String, DataBlock> map = hierarchicalStructure.getHierarchicalStructure();

			for (DataBlock db : map.values()) {
				List<Integer> index = db.getIndices();
				Matrix additionalGradient = new Matrix(getNumberOfParameters(), 1);	
				double inverseAdditionalLikelihoodTerm = 1d / getAdditionalLikelihoodTerm().get(index);			

				for (int i = 0; i < index.size() - 1; i++) {

					indexFirstObservation = index.get(i);
					setValuesInLikelihoodFunction(indexFirstObservation);
					likelihoodFirst =  getIndividualLikelihood().getValue();
					observedFirst =  getIndividualLikelihood().getYVector().getValueAt(0, 0);
					du_dbetaFirst =  getIndividualLikelihood().getGradient();

					for (int j = i + 1; j < index.size(); j++) {

						indexSecondObservation = index.get(j);
						setValuesInLikelihoodFunction(indexSecondObservation);
						likelihoodSecond =  getIndividualLikelihood().getValue();
						observedSecond =  getIndividualLikelihood().getYVector().getValueAt(0, 0);
						du_dbetaSecond =  getIndividualLikelihood().getGradient();

						sumObserved = observedFirst + observedSecond;

						multiplyingFactor = 1d;
						if (Math.abs(sumObserved - 1) < VERY_SMALL) {
							multiplyingFactor = -1d;
						}

						copulaExpression.setX(indexFirstObservation, indexSecondObservation);
						double copulaValue = copulaExpression.getValue();

						Matrix expansion1 = du_dbetaFirst.scalarMultiply((1 - likelihoodSecond) * -1d).add(
								du_dbetaSecond.scalarMultiply((1 - likelihoodFirst) * -1d)).scalarMultiply(copulaValue * multiplyingFactor * inverseAdditionalLikelihoodTerm);

						Matrix expansion2 = copulaExpression.getGradient().scalarMultiply((1 - likelihoodFirst) * 
								(1 - likelihoodSecond) * 
								multiplyingFactor * 
								inverseAdditionalLikelihoodTerm);

						tmp = expansion1.matrixStack(expansion2, true);

//						MatrixUtility.add(additionalGradient, tmp);
						additionalGradient = additionalGradient.add(tmp);

					}
				}

				additionalGradients.put(index, additionalGradient);

			}

			this.additionalGradients = additionalGradients;
			additionalGradientTermUptoDate = true;
		}
		
		return additionalGradients;
	}

	
	private Map<List<Integer>,SymmetricMatrix> getAdditionalHessians() {
		if (!additionalHessianTermUptoDate) {
			Map<List<Integer>, SymmetricMatrix> additionalHessians = new HashMap<List<Integer>, SymmetricMatrix>();
			
			Matrix additionalGradient;
				
			int indexFirstObservation;
			double likelihoodFirst;
			Matrix du_dbetaFirst;
			Matrix d2u_d2betaFirst;
			double observedFirst;

			int indexSecondObservation;
			double likelihoodSecond;
			Matrix du_dbetaSecond;
			Matrix d2u_d2betaSecond;
			double observedSecond;

			double sumObserved;
			double multiplyingFactor;

			Matrix tmp;

			Map<String, DataBlock> map = hierarchicalStructure.getHierarchicalStructure();
			for (DataBlock db : map.values()) {
				List<Integer> index = db.getIndices();
				Matrix additionalHessian = new Matrix(getNumberOfParameters(), getNumberOfParameters());

				double inverseAdditionalLikelihoodTerm = 1d / getAdditionalLikelihoodTerm().get(index);		
				additionalGradient = additionalGradients.get(index);

				additionalHessian.setSubMatrix(additionalGradient.multiply(additionalGradient.transpose()).scalarMultiply(-1d), 0, 0);	// first term corresponding to -1 * d1 ^ 2

				for (int i = 0; i < index.size() - 1; i++) {

					indexFirstObservation = index.get(i);
					setValuesInLikelihoodFunction(indexFirstObservation);
					likelihoodFirst =  getIndividualLikelihood().getValue();
					observedFirst =  getIndividualLikelihood().getYVector().getValueAt(0, 0);
					du_dbetaFirst =  getIndividualLikelihood().getGradient();
					d2u_d2betaFirst =  getIndividualLikelihood().getHessian();

					for (int j = i + 1; j < index.size(); j++) {

						indexSecondObservation = index.get(j);
						setValuesInLikelihoodFunction(indexSecondObservation);
						likelihoodSecond =  getIndividualLikelihood().getValue();
						observedSecond =  getIndividualLikelihood().getYVector().getValueAt(0, 0);
						du_dbetaSecond =  getIndividualLikelihood().getGradient();
						d2u_d2betaSecond =  getIndividualLikelihood().getHessian();

						sumObserved = observedFirst + observedSecond;

						multiplyingFactor = 1d;
						if (Math.abs(sumObserved - 1) < VERY_SMALL) {
							multiplyingFactor = -1d;
						}

						copulaExpression.setX(indexFirstObservation, indexSecondObservation);
						double copulaValue = copulaExpression.getValue();
						Matrix copulaGradient = copulaExpression.getGradient();
						Matrix copulaHessian = copulaExpression.getHessian();

						Matrix gradientMultipliedTemp = du_dbetaFirst.multiply(du_dbetaSecond.transpose());

						Matrix expansion11 = d2u_d2betaFirst.scalarMultiply((1 - likelihoodSecond) * -1d).add(
								d2u_d2betaSecond.scalarMultiply((1 - likelihoodFirst) * -1d)).add(							
										gradientMultipliedTemp.add(gradientMultipliedTemp.transpose())).scalarMultiply(
												copulaValue * multiplyingFactor * inverseAdditionalLikelihoodTerm);


						Matrix expansion12 = du_dbetaFirst.scalarMultiply((1 - likelihoodSecond) * -1d).add(
								du_dbetaSecond.scalarMultiply((1 - likelihoodFirst) * -1d)).multiply(copulaGradient.transpose()).scalarMultiply(multiplyingFactor * inverseAdditionalLikelihoodTerm);

						Matrix expansion22 = copulaHessian.scalarMultiply((1 - likelihoodFirst) * (1 - likelihoodSecond) * multiplyingFactor * inverseAdditionalLikelihoodTerm);

						tmp = expansion11.matrixStack(expansion12, false).matrixStack(expansion12.transpose().matrixStack(expansion22, false), true);

						additionalHessian = additionalHessian.add(tmp);
//						MatrixUtility.add(additionalHessian, tmp);

					}
				}

				additionalHessians.put(index, SymmetricMatrix.convertToSymmetricIfPossible(additionalHessian));
			}
			this.additionalHessians = additionalHessians;
			additionalHessianTermUptoDate = true;
		}
		return additionalHessians;
	}
	
	@Override
	public double getParameterValue(int index) {
		if (index < getIndividualLikelihood().getNumberOfParameters()) {
			return super.getParameterValue(index);
		} else if (index >= getIndividualLikelihood().getNumberOfParameters()) {
			index -= getIndividualLikelihood().getNumberOfParameters();
			return copulaExpression.getParameterValue(index);
		} else {
			return Double.NaN;
		}
	}

	@Override
	public int getNumberOfParameters() {
		return getIndividualLikelihood().getNumberOfParameters() + copulaExpression.getNumberOfParameters();
	}

	@Override
	public void setParameterValue(int index, double value) {
		if (index < getIndividualLikelihood().getNumberOfParameters()) {
			super.setParameterValue(index, value);
		} else if (index >= getIndividualLikelihood().getNumberOfParameters()) {
			index -= getIndividualLikelihood().getNumberOfParameters();
			copulaExpression.setParameterValue(index, value);
		} 
	}

	@Override
	public void optimizerDidThis(String actionString) {
		if (OptimizerListener.optimizationStarted.equals(actionString) || NewtonRaphsonOptimizer.InnerIterationStarted.equals(actionString)) {
			reset();
		}
	}

}

