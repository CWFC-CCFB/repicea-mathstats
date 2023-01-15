/*
 * This file is part of the repicea-simulation library.
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
package repicea.simulation.treelogger;

import java.awt.Container;

import javax.swing.JTextField;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import repicea.gui.REpiceaPanel;
import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;

/**
 * The TreeLogCategoryPanel class is the main class for any panel that allows the configuration
 * of a specific TreeLogCategory.
 * @author M. Fortin - August 2010
 */
@SuppressWarnings("serial")
public abstract class LogCategoryPanel<T extends LogCategory> extends REpiceaPanel implements CaretListener {

	/**
	 * This enum variable contains the common labels for a TreeLogCategoryPanel.
	 * @author Mathieu Fortin - April 2012
	 */
	public static enum MessageID implements TextableEnum {
		SmallEndDiameter("Minimum small-end diameter (cm)", "Diam\u00E8tre minimum au fin bout (cm)"),
		LongLength("Log length (m)", "Longueur de billon (m)"),
		LogGradeName("Log Category", "Cat\u00E9gorie de billon");

		
		MessageID(String englishText, String frenchText) {
			setText(englishText, frenchText);
		}
		

		@Override
		public void setText(String englishText, String frenchText) {
			REpiceaTranslator.setString(this, englishText, frenchText);
		}
		
		@Override
		public String toString() {return REpiceaTranslator.getString(this);}
		
	}

	
	protected T logCategory;

	protected JTextField nameTextField;

	/**
	 * General constructor for all the TreeLogCategoryPanel-derived classes.
	 */
	protected LogCategoryPanel(T logCategory) {
		super();
		this.logCategory = logCategory;
		nameTextField = new JTextField();
	}


	/**
	 * This method returns the name of the tree log category of this panel.
	 * @return a String
	 */
	public String getLogCategoryName() {return logCategory.getName();}
	
	/**
	 * This method returns the TreeLogCategory object this TreeLogCategoryPanel is based on. 
	 * @return a TreeLogCategory-derived instance
	 */
	protected T getTreeLogCategory() {return logCategory;}
	

	@Override
	public void caretUpdate(CaretEvent arg0) {
		if (arg0.getSource().equals(nameTextField)) {
			getTreeLogCategory().setName(nameTextField.getText());
			refreshUI();
		}
	}
	
	
	@SuppressWarnings("rawtypes")
	protected void refreshUI() {
		Container parent = getParent();
		while (!(parent instanceof TreeLoggerParametersDialog) && parent != null) {
			parent = parent.getParent();
		}
		if (parent instanceof TreeLoggerParametersDialog) {
			((TreeLoggerParametersDialog) parent).logCategoryList.repaint();
		}
	}
	
	@Override
	public void listenTo() {
		nameTextField.addCaretListener(this);
	}


	@Override
	public void doNotListenToAnymore() {
		nameTextField.removeCaretListener(this);
	}

	@Override
	public void refreshInterface() {
		nameTextField.setText(getLogCategoryName());
	}
	
}
