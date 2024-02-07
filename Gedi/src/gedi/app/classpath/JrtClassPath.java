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

import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.URL;
import java.nio.file.FileSystem;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.ProviderNotFoundException;

public class JrtClassPath implements ClassPath {

	private String root;
	private FileSystem fs;

	public static boolean hasJrt() {
		try {
			FileSystems.getFileSystem(URI.create("jrt:/"));
			return true;
		} catch (ProviderNotFoundException e) {
			return false;
		}
	}

	
	public JrtClassPath(String root) {
		try {
			this.fs = FileSystems.getFileSystem(URI.create("jrt:/"));
			this.root = root;
		} catch (ProviderNotFoundException e) {
			// Java 8!
		}
	}

	public boolean isResource(String resource) {
		if (fs==null) return false;
		return Files.isRegularFile(fs.getPath(root+resource));
	}

	public boolean isPackage(String packageName) {
		if (fs==null) return false;
		return Files.isDirectory(fs.getPath(root+packageName));
	}

	public String[] listPackages(String packageName) {
		if (fs==null) return new String[0];
		try {
			return Files.list(fs.getPath(root+packageName)).filter(p->Files.isDirectory(p)).map(JrtClassPath::cut).toArray(l->new String[l]);
		} catch (NoSuchFileException e) {
			return new String[0];
		} catch (IOException e) {
			throw new RuntimeException("Could not list jrt:",e);
		}
	}
	
	private static String cut(Path path) {
		return path.getFileName().toString();
	}

	public String[] listResources(String packageName) {
		if (fs==null) return new String[0];
		try {
			return Files.list(fs.getPath(root+packageName)).filter(p->Files.isRegularFile(p)).map(JrtClassPath::cut).toArray(l->new String[l]);
		} catch (NoSuchFileException e) {
			return new String[0];
		} catch (IOException e) {
			throw new RuntimeException("Could not list jrt:",e);
		}
	}
	
	@Override
	public URL getURL() {
		try {
			return URI.create("jrt:/").toURL();
		} catch (MalformedURLException e) {
			throw new RuntimeException();
		}
	}
	
	@Override
	public URL getResourceAsURL(String resource) throws MalformedURLException {
		while (resource.startsWith("/")) {
			resource = resource.substring(1);
		}
		return new URL(URI.create("jrt:/").toURL().toString()+resource);
	}

	public InputStream getResourceAsStream(String resource) {
		if (isResource(resource)) {
			try {
				return Files.newInputStream(fs.getPath(root+resource));
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		} else {
			return null;
		}
	}

	public String[] findResources(String rootPackageName, ResourceFilter resourceFilter) {
		return new ResourceFinder(this).findResources(rootPackageName, resourceFilter);
	}


	@Override
	public String toString() {
		return "jrt";
	}

	
}
