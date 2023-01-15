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

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

import javax.swing.AbstractButton;
import javax.swing.BorderFactory;
import javax.swing.JComponent;
import javax.swing.JOptionPane;
import javax.swing.KeyStroke;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import repicea.gui.CommonGuiUtility;
import repicea.gui.REpiceaAWTProperty;
import repicea.gui.REpiceaUIObject;
import repicea.gui.Resettable;
import repicea.gui.UIControlManager;
import repicea.gui.dnd.DnDPanel;
import repicea.gui.dnd.LocatedEvent;
import repicea.simulation.processsystem.UISetup.BasicMode;
import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;

@SuppressWarnings("serial")
public class SystemPanel extends DnDPanel<Processor> implements MouseListener, 
																ActionListener,
																Resettable, 
																ChangeListener {

	protected class SystemInternalPanel extends InternalPanel {
		
		@Override
		public void setSize(Dimension dim) {
			super.setSize(dim);
		}
		
		@Override
		public void paint(Graphics g) {
			super.paint(g);
			for (ValidProcessorLinkLine linkLine : linkLines) {
				linkLine.draw(g);
			}
			if (futureLink != null) {
				futureLink.draw(g);
			}
		}
		
	}
	
	
	private static enum MessageID implements TextableEnum {
		AboutToDeleteAProcessor("You are about to delete the selected processor or link. Are you sure you want to proceed?",
				"Vous allez effacer ce processus ou ce lien. Voulez-vous continuer?");
		
		MessageID(String englishText, String frenchText) {
			setText(englishText, frenchText);
		}

		@Override
		public void setText(String englishText, String frenchText) {
			REpiceaTranslator.setString(this, englishText, frenchText);
		}

		@Override
		public String toString() {
			return REpiceaTranslator.getString(this);
		}
	}

	protected final List<ValidProcessorLinkLine> linkLines;
	protected final List<ProcessorButton> processorButtons;
	@SuppressWarnings("rawtypes")
	private Enum currentMode;
	protected PreProcessorLinkLine futureLink;
	
	protected SystemPanel(SystemManager manager, SystemLayout systemLayout) {
		super(manager, Processor.class);
		setBorder(BorderFactory.createEtchedBorder());
		internalPanel.setLayout(systemLayout);
		
		linkLines = new CopyOnWriteArrayList<ValidProcessorLinkLine>();
		processorButtons = new CopyOnWriteArrayList<ProcessorButton>();
		
		ActionListener deleteListener = new ActionListener () {
			public void actionPerformed (ActionEvent actionEvent) {
				deleteAction();
			}
		};

		initUI();
		
		KeyStroke deleteStroke = KeyStroke.getKeyStroke(KeyEvent.VK_DELETE, 0);
		registerKeyboardAction(deleteListener, deleteStroke, JComponent.WHEN_IN_FOCUSED_WINDOW);
		KeyStroke backSpaceStroke = KeyStroke.getKeyStroke(KeyEvent.VK_BACK_SPACE, 0);
		registerKeyboardAction(deleteListener, backSpaceStroke, JComponent.WHEN_IN_FOCUSED_WINDOW);
	}
	
	
	protected void initUI() {
		internalPanel.removeAll();
		linkLines.clear();
		processorButtons.clear();
		futureLink = null;
		addManagerComponents();
		getListManager().checkForEndlessLoops();
		setMode(BasicMode.MoveProcessor);	// Default button
	}

	public SystemLayout getSystemLayout() {
		return (SystemLayout) internalPanel.getLayout();
	}
		
	@Override
	protected InternalPanel createInternalPanel() {return new SystemInternalPanel();}
	
	@Override
	protected void addManagerComponents() {
		for (REpiceaUIObject obj : getListManager().getList()) {
			Processor process = (Processor) obj;
			addProcessorButton(process.getUI(this));
			for (Processor subProcess : process.getSubProcessors()) {
				addLinkLine(new ProcessorLinkLine(this, process, subProcess));
			}
		}
	}

	
	protected boolean addProcessorButton(ProcessorButton processorButton) {
		boolean bool = false;
		if (!processorButtons.contains(processorButton)) {
			processorButtons.add(processorButton);
			add(processorButton);
			add(processorButton.label);
			processorButton.addMouseListener(this);
			bool = true;
		}
		return bool;
	}
	
	@Override
	public void acceptThisObject(Processor obj, LocatedEvent arg0) {
		Point dropLocation = getAbsoluteLocationFromDropEvent(arg0); // arg0.getLocation();
		ProcessorButton button = obj.getUI(this);
		Point originalLocation = new Point(dropLocation.x - button.getSize().width / 2, dropLocation.y - button.getSize().height / 2);
		obj.setOriginalLocation(originalLocation);
		if (addProcessorButton(button)) {		// light implementation to keep the synchronicity with the manager
			super.acceptThisObject(obj, arg0);
			SystemManagerDialog dlg = ((SystemManagerDialog) CommonGuiUtility.getParentComponent(this, SystemManagerDialog.class));
			if (dlg != null) {
				dlg.firePropertyChange(REpiceaAWTProperty.ActionPerformed, null, dlg);
			}
		}
	}

	protected void removeProcessorButton(ProcessorButton processorButton) {
		if (processorButtons.remove(processorButton)) {
			processorButton.removeMouseListener(this);
			processorButton.finalize();
			for (ValidProcessorLinkLine linkLine : linkLines) {
				if (linkLine.contains(processorButton)) {
					removeLinkLine(linkLine);
				}
			}
			getListManager().removeObject(processorButton.getOwner());
			SystemManagerDialog dlg = ((SystemManagerDialog) CommonGuiUtility.getParentComponent(this, SystemManagerDialog.class));
			if (dlg != null) {
				dlg.firePropertyChange(REpiceaAWTProperty.ActionPerformed, null, dlg);
			}
		}
	}

	@Override
	protected SystemManager getListManager() {return (SystemManager) super.getListManager();}
		
	@Override
	public void refreshInterface() {
		setMode(currentMode);
		for (ProcessorButton button : processorButtons) {
			button.setLabel();
		}
		internalPanel.getLayout().layoutContainer(this);
		internalPanel.revalidate();
		internalPanel.repaint();
	}
	
	protected void addLinkLine(ValidProcessorLinkLine linkLine) {
		if (!linkLines.contains(linkLine)) {
			linkLines.add(linkLine);
			linkLine.addMouseListener(this);
			SystemManagerDialog dlg = ((SystemManagerDialog) CommonGuiUtility.getParentComponent(this, SystemManagerDialog.class));
			if (dlg != null && linkLine.shouldChangeBeRecorded()) {
				dlg.firePropertyChange(REpiceaAWTProperty.ActionPerformed, null, dlg);
			}
		} else {
			linkLine.finalize();
		}
	}
	
	private void removeLinkLine(ValidProcessorLinkLine linkLine) {
		if (linkLines.remove(linkLine)) {
			linkLine.removeMouseListener(this);
			linkLine.finalize();
			SystemManagerDialog dlg = ((SystemManagerDialog) CommonGuiUtility.getParentComponent(this, SystemManagerDialog.class));
			if (dlg != null && linkLine.shouldChangeBeRecorded()) {
				dlg.firePropertyChange(REpiceaAWTProperty.ActionPerformed, null, dlg);
			}
		}
	}
	
	protected InternalPanel getInternalPanel() {
		return internalPanel;
	}
	
	private void deleteAction() {
		if (this.getListManager().getGUIPermission().isEnablingGranted()) {
			AbstractButton selectedButton = getSelectedFeature();
			if (selectedButton != null) {
				if (JOptionPane.showConfirmDialog(this, 
						MessageID.AboutToDeleteAProcessor.toString(), 
						UIControlManager.InformationMessageTitle.Warning.toString(), 
						JOptionPane.YES_NO_OPTION, 
						JOptionPane.WARNING_MESSAGE) == 0) {
					deleteFeature(selectedButton);
				}
			}
		}
	}

	protected void deleteFeature(AbstractButton button) {
		getViewport().setViewPositionVetoEnabled(true);	
		if (button instanceof ProcessorButton) {
			removeProcessorButton((ProcessorButton) button);
		} else if (button instanceof ValidProcessorLinkLine) {
			removeLinkLine((ValidProcessorLinkLine) button);
		}
		getListManager().checkForEndlessLoops();
		refreshInterface();
		sendVetoDisabledOnDispatchThread();
	}
	
	
	
	private AbstractButton getSelectedFeature() {
		AbstractButton selectedButton = null;
		List<AbstractButton> completeList = new ArrayList<AbstractButton>();
		completeList.addAll(processorButtons);
		completeList.addAll(linkLines);
		for (AbstractButton button : completeList) {
			if (button.isSelected()) {
				selectedButton = button;
				break;
			}
		}
		return selectedButton;
	}
	
	
	
	@Override
	public void mouseClicked(MouseEvent arg0) {
		List<AbstractButton> completeList = new ArrayList<AbstractButton>();
		completeList.addAll(processorButtons);
		completeList.addAll(linkLines);
		for (AbstractButton button : completeList) {
			button.setSelected(button.equals(arg0.getSource()));
		}
		internalPanel.revalidate();
		internalPanel.repaint();
	}

	@Override
	public void mouseEntered(MouseEvent arg0) {}

	@Override
	public void mouseExited(MouseEvent arg0) {}

	@Override
	public void mousePressed(MouseEvent arg0) {}

	@Override
	public void mouseReleased(MouseEvent arg0) {}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		if (arg0.getSource() instanceof ToolButton) {
			ToolButton button = (ToolButton) arg0.getSource();
			setMode(button.getMode());
		}
	}

	@SuppressWarnings("rawtypes")
	void setMode(Enum mode) {
		if (mode != null) {
			currentMode = mode; 
			for (ProcessorButton button : processorButtons) {
				button.setDragMode(currentMode);
			}
		}
	}

	Enum getMode() {return currentMode;}
	
	protected void registerLinkBeingCreated(PreProcessorLinkLine futureLink) {
		if (this.futureLink != null) {
			this.futureLink.finalize();
		} 
		this.futureLink = futureLink;
	}

	protected void doLinkAction() {
		if (!linkLines.contains(futureLink)) {
			addLinkLine(futureLink.convertIntoProcessorLinkLine());
		}
		registerLinkBeingCreated(null);
		getListManager().checkForEndlessLoops();
		validate();
	}

	@Override
	public void reset() {initUI();}

	@Override
	public void stateChanged(ChangeEvent e) {
		if (e.getSource().equals(getListManager().guiInterface.zoomSlider)) {
			int newValue = getListManager().guiInterface.zoomSlider.getValue();
			getSystemLayout().setCurrentZoom(newValue);
			refreshInterface();
		}
	}


}
