/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2021 Mathieu Fortin for Rouge-Epicea
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

import java.awt.Point;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.SwingUtilities;

import org.junit.Assert;
import org.junit.Test;

import repicea.gui.REpiceaAWTProperty;
import repicea.gui.REpiceaGUITestRobot;
import repicea.gui.UIControlManager;
import repicea.gui.UIControlManager.CommonControlID;
import repicea.gui.dnd.LocatedEvent;
import repicea.simulation.processsystem.SystemManagerDialog.MessageID;
import repicea.util.ObjectUtility;

public class ProcessSystemTest {

//	private final static int WAIT_TIME = 500;

	private static enum FakeEnum {class1, class2}
	
	/**
	 * This test loads a file and checks if the process system can be read.
	 */
	@Test
	public void loadSystemTest() {
		String filename = ObjectUtility.getPackagePath(getClass()) + "systemTest1.prl";
		SystemManager manager = new SystemManager();
		try {
			manager.load(filename);
			System.out.println("Process system systemTest1.prl read!");
		} catch (IOException e) {
			Assert.fail();
		}
	}

	
	/**
	 * This test if the amount in input is equal to the sum of the amounts in output.
	 */
	@SuppressWarnings("rawtypes")
	@Test
	public void matterBalanceTest() {
		String filename = ObjectUtility.getPackagePath(getClass()) + "systemTest1.prl";
		SystemManager manager = new SystemManager();
		try {
			manager.load(filename);
			System.out.println("Process system systemTest1.prl read!");
		} catch (IOException e) {
			Assert.fail();
		}

		List<Processor> primaryProcessors = manager.getPrimaryProcessors();
		ProcessUnit<FakeEnum> inputUnit = new ProcessUnit<FakeEnum>();
		inputUnit.getAmountMap().put(FakeEnum.class1, 100d);
		inputUnit.getAmountMap().put(FakeEnum.class2, 100d);
		List<ProcessUnit> processUnits = new ArrayList<ProcessUnit>();
		processUnits.add(inputUnit);
		Collection<ProcessUnit> outputProcessUnits = primaryProcessors.get(0).doProcess(processUnits);
		double sumClass1 = 0d;
		double sumClass2 = 0d;
		for (ProcessUnit processUnit : outputProcessUnits) {
			sumClass1 += (Double) processUnit.getAmountMap().get(FakeEnum.class1);
			sumClass2 += (Double) processUnit.getAmountMap().get(FakeEnum.class2);
		}
		Assert.assertEquals(100d, sumClass1, 1E-10);
		Assert.assertEquals(100d, sumClass2, 1E-10);
	}
	

	@Test
	public void creatingProcessorsLinksAndChangingValuesOfSliderAsASingleAction() throws InterruptedException {
		SystemManager man = new SystemManager();
		SystemManagerDialog dlg = man.getUI(null);
		
		REpiceaGUITestRobot robot = new REpiceaGUITestRobot();
		Thread t = robot.showWindow(man);
		
		SystemPanel sysPane = dlg.systemPanel;
		Point loc = sysPane.getLocation();
		LocatedEvent evt = new LocatedEvent(sysPane, new Point(loc.x + 50, loc.y + 50));
			sysPane.acceptThisObject(new Processor("1"), evt);
		REpiceaGUITestRobot.letDispatchThreadProcess();
		Assert.assertEquals("Testing if the processor has been recorded in the manager", 1, man.getList().size());
		
		evt = new LocatedEvent(sysPane, new Point(loc.x + 125, loc.y + 50));
		sysPane.acceptThisObject(new Processor("2"), evt);
		REpiceaGUITestRobot.letDispatchThreadProcess();
		Assert.assertEquals("Testing if the processor has been recorded in the manager", 2, man.getList().size());
		
		sysPane.addLinkLine(new ProcessorLinkLine(sysPane, man.getList().get(0), man.getList().get(1)));
		REpiceaGUITestRobot.letDispatchThreadProcess();
		Assert.assertEquals("Testing if the link has been recorded in the manager", 1, sysPane.linkLines.size());
		Assert.assertEquals("Testing if the link has been recorded in father processor", 1, man.getList().get(0).getSubProcessors().size());

		dlg.undo.doClick();
		REpiceaGUITestRobot.letDispatchThreadProcess();
		Assert.assertEquals("Testing if the link has been undone in the manager", 0, sysPane.linkLines.size());
		Assert.assertEquals("Testing if the link has been undone in father processor", 0, man.getList().get(0).getSubProcessors().size());
		
		dlg.redo.doClick();
		REpiceaGUITestRobot.letDispatchThreadProcess();
		Assert.assertEquals("Testing if the link has been redone in the manager", 1, sysPane.linkLines.size());
		Assert.assertEquals("Testing if the link has been redone in father processor", 1, man.getList().get(0).getSubProcessors().size());

		final ProcessorLinkLine link = (ProcessorLinkLine) sysPane.linkLines.get(0);
		int intake = man.getList().get(0).getSubProcessorIntakes().get(man.getList().get(1));
		Assert.assertEquals("Testing if initial flow is set to 0", 0, intake);
		ProcessorLinkLineSlider linkDlg = link.getUI(dlg);
		
		Thread t2 = robot.showWindow(link);
		
		REpiceaGUITestRobot.letDispatchThreadProcess();
		linkDlg.slider.setValue(25);
		REpiceaGUITestRobot.letDispatchThreadProcess();
		linkDlg.slider.setValue(55);
		REpiceaGUITestRobot.letDispatchThreadProcess();
		
		linkDlg.setVisible(false);
		linkDlg.windowClosing(new WindowEvent(linkDlg, WindowEvent.WINDOW_CLOSING));
		REpiceaGUITestRobot.letDispatchThreadProcess();

		intake = man.getList().get(0).getSubProcessorIntakes().get(man.getList().get(1));
		Assert.assertEquals("Testing if flow is now set to 55", 55, intake);
		
		dlg.undo.doClick();
		REpiceaGUITestRobot.letDispatchThreadProcess();
		intake = man.getList().get(0).getSubProcessorIntakes().get(man.getList().get(1));
		Assert.assertEquals("Testing if intake has been properly undone", 0, intake);

		dlg.redo.doClick();
		REpiceaGUITestRobot.letDispatchThreadProcess();
		intake = man.getList().get(0).getSubProcessorIntakes().get(man.getList().get(1));
		Assert.assertEquals("Testing if intake has been properly redone", 55, intake);

		t2.join();
		dlg.setVisible(false);
		dlg.dispose();
		t.join();
	}

	
	@Test
	public void testSVNExportHappyPath() throws Exception {
		SystemManager man = new SystemManager();
		SystemManagerDialog dlg = man.getUI(null);
		
		REpiceaGUITestRobot robot = new REpiceaGUITestRobot();
		Thread t = robot.showWindow(man);
		
		robot.clickThisButton(UIControlManager.CommonMenuTitle.File.name()); 
		
		Runnable toRun = new Runnable() {
			public void run() {
				try {
					robot.clickThisButton(MessageID.ExportToSVG.name());
				} catch (InvocationTargetException | InterruptedException e) {} 
			}
		};
		robot.startGUI(toRun, JDialog.class);

		String filename = System.getProperty("java.io.tmpdir") + "exportFluxConfiguration.svg";
		File f = new File(filename);
		if (f.exists())
			f.delete();
		robot.fillThisTextField("Filename", filename);
		
		robot.clickThisButton(CommonControlID.Save.name(), REpiceaAWTProperty.SVGFileSaved);
		
		Assert.assertTrue("Testing if svg export file exists", f.exists());

		robot.clickThisButton(UIControlManager.CommonMenuTitle.File.name()); 
		robot.clickThisButton(UIControlManager.CommonControlID.Close.name());
	}


	
}
