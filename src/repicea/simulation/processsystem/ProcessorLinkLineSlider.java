/*
 * This file is part of the repicea-simulation library.
 *
 * Copyright (C) 2009-2014 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation.processsystem;

import java.awt.FlowLayout;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.util.Map;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import repicea.gui.CommonGuiUtility;
import repicea.gui.REpiceaAWTProperty;
import repicea.gui.REpiceaDialog;
import repicea.gui.UIControlManager;
import repicea.simulation.processsystem.SystemManagerDialog.MessageID;

@SuppressWarnings("serial")
public class ProcessorLinkLineSlider extends REpiceaDialog implements ChangeListener, ActionListener, WindowListener {

	private final ProcessorLinkLine linkLine;
	protected final JSlider slider;
	protected final JButton setToAppropriateValue;
	private transient int startingValue;
	
	protected ProcessorLinkLineSlider(Window window, ProcessorLinkLine linkLine) {
		super(window);
		this.linkLine = linkLine;
		slider = new JSlider();
		slider.setMaximum(100);
		slider.setMinimum(0);
		slider.setMajorTickSpacing(20);
		slider.setPaintTicks(true);
		slider.setPaintLabels(true);
		setToAppropriateValue = UIControlManager.createButtonWithRedCircleIcon();
		initUI();
		refreshInterface();
		pack();
		setResizable(false);
	}
		
	@Override
	protected void initUI() {
		setTitle(MessageID.SliderTitle.toString());
		getContentPane().setLayout(new FlowLayout(FlowLayout.LEADING));
		getContentPane().add(Box.createHorizontalStrut(5));
		getContentPane().add(slider);
		getContentPane().add(Box.createHorizontalStrut(10));
		getContentPane().add(setToAppropriateValue);
		getContentPane().add(Box.createHorizontalStrut(5));
	}
	
	
	@Override
	public void listenTo() {
		slider.addChangeListener(this);
		setToAppropriateValue.addActionListener(this);
		addWindowListener(this);
	}

	@Override
	public void doNotListenToAnymore() {
		slider.removeChangeListener(this);
		setToAppropriateValue.removeActionListener(this);
		removeWindowListener(this);
	}

	@Override
	public void refreshInterface() {
		Processor fatherProcessor = linkLine.getFatherAnchor().getOwner();
		Processor sonProcessor = linkLine.getSonAnchor().getOwner();
		int value = fatherProcessor.getSubProcessorIntakes().get(sonProcessor);
		slider.setValue(value);
		super.refreshInterface();
	}
	
	
	@Override
	public void stateChanged(ChangeEvent arg0) {
		if (arg0.getSource().equals(slider)) {
			int value = slider.getValue();
			Processor fatherProcessor = linkLine.getFatherAnchor().getOwner();
			Processor sonProcessor = linkLine.getSonAnchor().getOwner();
			fatherProcessor.getSubProcessorIntakes().put(sonProcessor, value);
			linkLine.setLabel();
//			SystemManagerDialog dlg = (SystemManagerDialog) CommonGuiUtility.getParentComponent(this, SystemManagerDialog.class);
//			if (dlg != null) {
//				dlg.firePropertyChange(REpiceaAWTProperty.ActionPerformed, null, dlg);
//			}
		}
	}

	private int getSumOfCurrentFlux() {
		Processor fatherProcessor = linkLine.getFatherAnchor().getOwner();
		Processor sonProcessor = linkLine.getSonAnchor().getOwner();
		int sum = 0;
		Map<Processor, Integer> processorMap = fatherProcessor.getSubProcessorIntakes();
		for (Processor subProcessor : processorMap.keySet()) {
			if (!subProcessor.equals(sonProcessor)) {
				sum += processorMap.get(subProcessor);
			}
		}
		return sum;
	}
	
	
	@Override
	public void actionPerformed(ActionEvent arg0) {
		if (arg0.getSource().equals(setToAppropriateValue)) {
			Processor fatherProcessor = linkLine.getFatherAnchor().getOwner();
			Processor sonProcessor = linkLine.getSonAnchor().getOwner();
//			int originalValue = fatherProcessor.getSubProcessorIntakes().get(sonProcessor);
			int newValue = 100 - getSumOfCurrentFlux();
			fatherProcessor.getSubProcessorIntakes().put(sonProcessor, newValue);
			linkLine.setLabel();
			refreshInterface();
//			SystemManagerDialog dlg = (SystemManagerDialog) CommonGuiUtility.getParentComponent(this, SystemManagerDialog.class);
//			if (dlg != null) {
//				dlg.firePropertyChange(REpiceaAWTProperty.ActionPerformed, originalValue, newValue);
//			}
		}
	}

	@Override
	public void windowOpened(WindowEvent e) {
		startingValue = getCurrentValue();
	}
	
	
	private int getCurrentValue() {
		Processor fatherProcessor = linkLine.getFatherAnchor().getOwner();
		Processor sonProcessor = linkLine.getSonAnchor().getOwner();
		return fatherProcessor.getSubProcessorIntakes().get(sonProcessor);
	}

	@Override
	public void windowClosing(WindowEvent e) {
		int currentValue = getCurrentValue();
		if (currentValue != startingValue) {
			SystemManagerDialog dlg = (SystemManagerDialog) CommonGuiUtility.getParentComponent(this, SystemManagerDialog.class);
			if (dlg != null) {
				dlg.firePropertyChange(REpiceaAWTProperty.ActionPerformed, startingValue, currentValue);
			}
		}
	}

	@Override
	public void windowClosed(WindowEvent e) {}

	@Override
	public void windowIconified(WindowEvent e) {}

	@Override
	public void windowDeiconified(WindowEvent e) {}

	@Override
	public void windowActivated(WindowEvent e) {}

	@Override
	public void windowDeactivated(WindowEvent e) {}


	
	
	
}
