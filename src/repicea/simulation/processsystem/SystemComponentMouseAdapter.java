/*
 * This file is part of the repicea-simulation library.
 *
 * Copyright (C) 2009-2015 Mathieu Fortin for Rouge-Epicea
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

import java.awt.Component;
import java.awt.Container;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.security.InvalidParameterException;

import javax.swing.AbstractButton;
import javax.swing.JComboBox;
import javax.swing.JSlider;
import javax.swing.SwingUtilities;
import javax.swing.text.JTextComponent;

import repicea.gui.CommonGuiUtility;
import repicea.gui.REpiceaShowableUIWithParent;
import repicea.gui.Refreshable;
import repicea.gui.permissions.REpiceaGUIPermissionProvider;

class SystemComponentMouseAdapter extends MouseAdapter {

	private final REpiceaGUIPermissionProvider readWriteProvider;
	
	
	SystemComponentMouseAdapter(REpiceaGUIPermissionProvider readWriteProvider) {
		if (readWriteProvider instanceof REpiceaShowableUIWithParent && readWriteProvider instanceof Component) {
			this.readWriteProvider = readWriteProvider;
		} else {
			throw new InvalidParameterException("The component must implement the ShowableObjectWithParent interface");
		} 
	}
	
	
	
	
	@Override
	public void mouseClicked(MouseEvent arg0) {
		REpiceaShowableUIWithParent showable = (REpiceaShowableUIWithParent) readWriteProvider;
		if (SwingUtilities.isLeftMouseButton(arg0) && arg0.getClickCount() >= 2) {
			SystemManagerDialog systemManagerDialog = (SystemManagerDialog) CommonGuiUtility.getParentComponent((Component) readWriteProvider, SystemManagerDialog.class);
			boolean overallEnabling = systemManagerDialog.getCaller().getGUIPermission().isEnablingGranted();
			boolean isEnablingGranted = readWriteProvider.getGUIPermission().isEnablingGranted() && overallEnabling;
			Container internalDlg = (Container) showable.getUI(systemManagerDialog);
			if (internalDlg instanceof Refreshable) {
				((Refreshable) internalDlg).refreshInterface();
			}
			if (!isEnablingGranted) {
				CommonGuiUtility.enableThoseComponents(internalDlg, JTextComponent.class, isEnablingGranted);
				CommonGuiUtility.enableThoseComponents(internalDlg, JSlider.class, isEnablingGranted);
				CommonGuiUtility.enableThoseComponents(internalDlg, AbstractButton.class, isEnablingGranted);
				CommonGuiUtility.enableThoseComponents(internalDlg, JComboBox.class, isEnablingGranted);
			}
			showable.showUI(systemManagerDialog);
		}
	}

}
