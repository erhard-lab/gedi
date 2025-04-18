/*
 * Copyright 2007 Google Inc.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy of
 * the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations under
 * the License.
 */
package gedi.app.classpath;

import java.io.File;
import java.util.ArrayList;

public class ClassPathFactory {

	public static final String JAVA_CLASS_PATH = "java.class.path";

	private static ClassPathFactory instance = null;
	public static ClassPathFactory getInstance() { 
		if (instance==null) instance = new ClassPathFactory();
		return instance;
	}


	private ClassPathFactory() {}

	public String getJVMClasspath() {
		return System.getProperty(JAVA_CLASS_PATH);
	}

	public String[] parseClasspath(String classpath) {
		return classpath.split(File.pathSeparator);
	}

	public CompositeClassPath createFromJVM() {
		return createFromPath(getJVMClasspath());
	}

	public CompositeClassPath createFromPath(String classpath) {
		return createFromPaths(parseClasspath(classpath));
	}

	public CompositeClassPath createFromPaths(String... paths) {
		ArrayList<ClassPath> classPaths = new ArrayList<ClassPath>();
		for (String path : paths) {
			File file = new File(path);
			if (file.isFile()) {
				classPaths.add(new JARClassPath(file));
			} else if (file.isDirectory()) {
				classPaths.add(new DirectoryClassPath(file));
			} else if (path.equals("jrt") && JrtClassPath.hasJrt()) {
				classPaths.add(new JrtClassPath("/modules/java.base/"));
				classPaths.add(new JrtClassPath("/modules/java.desktop/"));
			} else {
				//ignore since that is what JVM does.
			}
		}
		ClassPath[] array = new ClassPath[classPaths.size()];
		array = (ClassPath[]) classPaths.toArray(array);
		return new CompositeClassPath(array);
	}

}
