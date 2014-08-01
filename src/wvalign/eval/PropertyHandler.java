package wvalign.eval;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

public class PropertyHandler {

	private String[] requiredProps;

	public PropertyHandler(String[] requires) {
		requiredProps = requires;
	}

	public boolean checkProps(Properties properties) {
		for (String reqProp : requiredProps) {
			if (!properties.containsKey(reqProp)) {
				System.out.println("Error! Missing property: " + reqProp);
				return false;
			}
		}
		return true;
	}

	public Properties readPropertiesFile(String filename) {
		Properties props = new Properties();
		try {
			props.load(new FileInputStream(filename));
		} catch (IOException ex) {
			System.out.println("Could not open properties file " + filename);
		}
		return props;
	}
}
