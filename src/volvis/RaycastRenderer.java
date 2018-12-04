/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        //tFunc.setTestFunc();
        
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     
    int count = 0;
    short getVoxel(double[] coord) {
       
        if (coord[0] < 0 || Math.ceil(coord[0]) >= volume.getDimX() || coord[1] < 0 || Math.ceil(coord[1]) >= volume.getDimY()
                || coord[2] < 0 || Math.ceil(coord[2]) >= volume.getDimZ()) {
            return 0;
        }
        
        int x = (int) Math.ceil(coord[0]);
        int y = (int) Math.ceil(coord[1]);
        int z = (int) Math.ceil(coord[2]);

        return volume.getVoxel(x, y, z);
    }
    
    short getTrilinearVoxel(double[] coord) {
        if (coord[0] < 0 || Math.ceil(coord[0]) >= volume.getDimX() || coord[1] < 0 || Math.ceil(coord[1]) >= volume.getDimY()
                || coord[2] < 0 || Math.ceil(coord[2]) >= volume.getDimZ()) {
            return 0;
        }
        
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        
        double alpha, alpha_expression;
        double beta, beta_expression;
        double gamma, gamma_expression;

        alpha_expression = (x - Math.floor(x))/(Math.ceil(x) - Math.floor(x));
        beta_expression = (y - Math.floor(y))/(Math.ceil(y) - Math.floor(y));
        gamma_expression = (z - Math.floor(z))/(Math.ceil(z) - Math.floor(z));

        alpha = (!Double.isNaN((alpha_expression))) ? alpha_expression : 0.0;
        beta = (!Double.isNaN(beta_expression)) ? beta_expression : 0.0;
        gamma = (!Double.isNaN(gamma_expression)) ? gamma_expression : 0.0;
                
        double x0 = volume.getVoxel((int) Math.floor(x), (int) Math.floor(y), (int) Math.floor(z));
        double x1 = volume.getVoxel((int) Math.ceil(x), (int) Math.floor(y), (int) Math.floor(z));
        double x2 = volume.getVoxel((int) Math.floor(x), (int) Math.floor(y), (int) Math.ceil(z));
        double x3 = volume.getVoxel((int) Math.ceil(x), (int) Math.floor(y), (int) Math.ceil(z));
        double x4 = volume.getVoxel((int) Math.floor(x), (int) Math.ceil(y), (int) Math.floor(z));
        double x5 = volume.getVoxel((int) Math.ceil(x), (int) Math.ceil(y), (int) Math.floor(z));
        double x6 = volume.getVoxel((int) Math.floor(x), (int) Math.ceil(y), (int) Math.ceil(z));
        double x7 = volume.getVoxel((int) Math.ceil(x), (int) Math.ceil(y), (int) Math.ceil(z));
   
        
        double voxelVal = (1 - alpha)*(1 - beta)*(1 - gamma)*x0 + alpha*(1 - beta)*(1 - gamma)*x1
                + (1 - alpha)*beta*(1 - gamma)*x2 + alpha*beta*(1 - gamma)*x3
                + (1 - alpha)*(1 - beta)*gamma*x4 + alpha*(1 - beta)*gamma*x5
                + (1 - alpha)*beta*gamma*x6 + alpha*beta*gamma*x7;
        
        return (short)voxelVal;

    }

    void slicer(double[] viewMatrix) {
         
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                int val = getVoxel(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }
    
    void mip(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                int maxVal = 0;
                for (int k = 0; k < volume.getDimZ() - 1; k+=1) {
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + viewVec[0] * (k - volumeCenter[0]) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + viewVec[1] * (k - volumeCenter[1]) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + viewVec[2] * (k - volumeCenter[2]) + volumeCenter[2];
                    // Get trilinearVoxel(pixelCoord) instead of getVoxel(pixelCoord) after implememtning trilinearlization
                    if (getTrilinearVoxel(pixelCoord) > maxVal) {
                        maxVal = getTrilinearVoxel(pixelCoord);
                    }
                }
                int val = maxVal;
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }
        
    }
    
    void composite(double[] viewMatrix) {
         
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor;
        ArrayList<TFColor> colors;
        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {

                colors = new ArrayList<TFColor>();
                for (int k = 0; k < volume.getDimZ() - 1; k+=15) {
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + viewVec[0] * (k - volumeCenter[0]) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + viewVec[1] * (k - volumeCenter[1]) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + viewVec[2] * (k - volumeCenter[2]) + volumeCenter[2];
                    // Get trilinearVoxel(pixelCoord) instead of getVoxel(pixelCoord) after implememtning trilinearlization
                    TFColor color;
                            
                    color = getColorInterpolated(pixelCoord);
                    colors.add(color);
                }
                
                double red = getAverageColor('r', 0, colors);
                double blue = getAverageColor('b', 0, colors);
                double green = getAverageColor('g', 0, colors);
                   
//                System.out.printf("%f %f %f\n", red, green, blue);
                
                double[] colorList = new double[]{red, green, blue};
                voxelColor = new TFColor(colorList[0], colorList[1], colorList[2], 1);
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }
    }
    
    public boolean alreadyInteger(double[] coord) {
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];

        double[] xDim = {Math.floor(x), Math.ceil(x)};
        double[] yDim = {Math.floor(y), Math.ceil(y)};
        double[] zDim = {Math.floor(z), Math.ceil(z)};

        if (xDim[0] == xDim[1] && yDim[0] == yDim[1] && zDim[0] == zDim[1]) {
            return true;
        }
        return false;
    }
    
    public TFColor getColorInterpolated(double[] coord) {
        if (alreadyInteger(coord)) {
            return tFunc.getColor(getVoxel(coord));
        }
        
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];

        double[] xDim = {Math.floor(x), Math.ceil(x)};
        double[] yDim = {Math.floor(y), Math.ceil(y)};
        double[] zDim = {Math.floor(z), Math.ceil(z)};
       
        TFColor[] colors = new TFColor[8];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    for (int index = 0; index < colors.length; index++) {
                        double[] voxelCoord = {xDim[i], yDim[j], zDim[k]};
                        colors[index] = tFunc.getColor(getVoxel(voxelCoord));
                    }
                }
            }
        }
        
        double[][] colorsNeighboringVoxels = getColorValues(colors);

        double[] red = colorsNeighboringVoxels[0];
        double[] green = colorsNeighboringVoxels[1];
        double[] blue = colorsNeighboringVoxels[2];
        double[] alpha = colorsNeighboringVoxels[3];

        TFColor result = new TFColor(
                getColorInterpolated(red, coord),
                getColorInterpolated(green, coord),
                getColorInterpolated(blue, coord),
                getColorInterpolated(alpha, coord));
        return result;
    }
    
    public double[][] getColorValues(TFColor[] colors) {
        double[] red = new double[colors.length];
        double[] green = new double[colors.length];
        double[] blue = new double[colors.length];  
        double[] alpha = new double[colors.length];

        for (int i = 0; i < colors.length; i++) {
            red[i] = colors[i].r;
            green[i] = colors[i].g;
            blue[i] = colors[i].b;
            alpha[i] = colors[i].a;
        }
        return new double[][]{red, green, blue, alpha};
    }
    
    public double getColorInterpolated(double[] colors, double[] coord) {
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];

        double[] xDim = {Math.floor(x), Math.ceil(x)};
        double[] yDim = {Math.floor(y), Math.ceil(y)};
        double[] zDim = {Math.floor(z), Math.ceil(z)};
        
        double xd = x - xDim[0] / (xDim[1] - xDim[0]);
        double yd = y - yDim[0] / (yDim[1] - yDim[0]);
        double zd = z - zDim[0] / (zDim[1] - zDim[0]);
        
        double C00 = colors[0] * (1 - xd) + colors[4] * xd;
        double C01 = colors[2] * (1 - xd) + colors[6] * xd;
        double C10 = colors[1] * (1 - xd) + colors[5] * xd;
        double C11 = colors[3] * (1 - xd) + colors[7] * xd;

        double C0 = C00 * (1 - yd) + C10 * yd;
        double C1 = C01 * (1 - yd) + C11 * yd;

        double C = C0 * (1 - zd) + C1 * zd;
        return C;
    }
    
    double getAverageColor(char color, int voxel, ArrayList<TFColor> colors) {
        if (voxel >= colors.size()) {
            return 0;
        }
        if (color == 'b') {
            return colors.get(voxel).a * colors.get(voxel).b + (1 - colors.get(voxel).a) * getAverageColor(color, voxel + 1, colors);
        }
        if (color == 'r') {
            return colors.get(voxel).a * colors.get(voxel).r + (1 - colors.get(voxel).a) * getAverageColor(color, voxel + 1, colors);
        }
        else {
            return colors.get(voxel).a * colors.get(voxel).g + (1 - colors.get(voxel).a) * getAverageColor(color, voxel + 1, colors);
        }
        
    }
    
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
//        
//        slicer(viewMatrix);
//        mip(viewMatrix);
        composite(viewMatrix);
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
