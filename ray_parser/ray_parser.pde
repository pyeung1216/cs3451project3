///////////////////////////////////////////////////////////////////////
//
// Command Line Interface (CLI) Parser  
//
///////////////////////////////////////////////////////////////////////
String gCurrentFile = new String("rect_test.cli"); // A global variable for holding current active file name.
float angle, k, bgR, bgG, bgB;

ArrayList<Light> lights;
ArrayList<Shape> shapes;

Surface surface;
Triangle triangle;

///////////////////////////////////////////////////////////////////////
//
// Classes
//
///////////////////////////////////////////////////////////////////////

class Vertex {
  float[] v;
  
  Vertex()  {
    v = new float[3];
  }
  
  Vertex(float x, float y, float z) {
    v = new float[3];
    v[0] = x;
    v[1] = y;
    v[2] = z;
  }
  
  Vertex add(Vertex add)  {
    Vertex added = new Vertex();
    added.v[0] = v[0] + add.v[0];
    added.v[1] = v[1] + add.v[1];
    added.v[2] = v[2] + add.v[2];
    return added;
  }
  
  Vertex subtract(Vertex sub)  {
    Vertex subtracted = new Vertex();
    subtracted.v[0] = v[0] - sub.v[0];
    subtracted.v[1] = v[1] - sub.v[1];
    subtracted.v[2] = v[2] - sub.v[2];
    return subtracted;
  }
  
  float dotProduct(Vertex dot)  {
    return v[0] * dot.v[0] + v[1] * dot.v[1] + v[2] * dot.v[2];
  }
  
  Vertex crossProduct(Vertex cross)  {
    Vertex crossed = new Vertex();
    crossed.v[0] = v[1] * cross.v[2] - v[2] * cross.v[1];
    crossed.v[1] = v[2] * cross.v[0] - v[0] * cross.v[2];
    crossed.v[2] = v[0] * cross.v[1] - v[1] * cross.v[0];
    return crossed;
  }
  
  float euclideanDistance(Vertex v2)  {
    float euclidDist = sqrt((pow((v[0] - v2.v[0]), 2)) +
                            (pow((v[1] - v2.v[1]), 2)) + 
                            (pow((v[2] - v2.v[2]), 2)));
    return euclidDist;
  }
    
  float getX(){
    return v[0];
  }
  
  float getY(){
    return v[1];
  }
  
  float getZ(){
    return v[2];
  }
}

class Ray  {
  Vertex origin;
  Vertex vector;
  
  Ray()  {
    origin = new Vertex();
    vector = new Vertex();
  }
  
  void origin(float x, float y, float z)  {
    origin.v[0] = x;
    origin.v[1] = y;
    origin.v[2] = z;
  }
  
  void vector(float x, float y, float z)  {
    vector.v[0] = x;
    vector.v[1] = y;
    vector.v[2] = z;
  }
  
  Vertex getOrigin()  {
    return origin;
  }
  float getOriginX(){
    return origin.v[0];
  }
  float getOriginY(){
    return origin.v[1];
  }
  float getOriginZ(){
    return origin.v[2];
  }
  
  Vertex getVector()  {
    return vector;
  }
  float getVectorX(){
    return vector.v[0];
  }
  float getVectorY(){
    return vector.v[1];
  }
  float getVectorZ(){
    return vector.v[2];
  }
}

abstract class Shape  {
  Surface shapeDetails;
  
  void surface(Surface s)  {
    shapeDetails = s;
  }
  
  abstract float intersects(Ray r);
  abstract Vertex getIntersect(Ray r, float t);
  abstract Vertex getIntersectNorm(Ray r, float t);
}

class Triangle extends Shape  {
  ArrayList<Vertex> triangle;
  
  Triangle()  {
    triangle = new ArrayList<Vertex>();
  }
  
  void vertex(float x, float y, float z)  {    
    Vertex v = new Vertex(x, y, z);
    triangle.add(v);
  }   
  
  float intersects(Ray r)  {
    Vertex edge1 = triangle.get(1).subtract(triangle.get(0));
    Vertex edge2 = triangle.get(2).subtract(triangle.get(1));    
    Vertex edge3 = triangle.get(0).subtract(triangle.get(2));
    Vertex normal = edge1.crossProduct(edge2);
    
    float normMag = sqrt(pow(normal.v[0], 2) + pow(normal.v[1], 2) + pow(normal.v[2], 2));
    normal.v[0] = normal.v[0] / normMag;
    normal.v[1] = normal.v[1] / normMag;
    normal.v[2] = normal.v[2] / normMag;
    
    float dVal = -1 * (normal.v[0] * triangle.get(0).v[0] + normal.v[1] * triangle.get(0).v[1] + normal.v[2] * triangle.get(0).v[2]);
    
    if(normal.dotProduct(r.getVector()) == 0)  {
      return -1;
    }
    
    float tVal = -1 * (normal.v[0] * r.getOrigin().v[0] + normal.v[1] * r.getOrigin().v[1] + normal.v[2] * r.getOrigin().v[2] + dVal) / (normal.dotProduct(r.getVector()));
    
    float xIntersect = r.getOrigin().v[0] + tVal * r.getVector().v[0];
    float yIntersect = r.getOrigin().v[1] + tVal * r.getVector().v[1];
    float zIntersect = r.getOrigin().v[2] + tVal * r.getVector().v[2];
    Vertex intersect = new Vertex(xIntersect, yIntersect, zIntersect);
    
    Vertex P0Intersect = intersect.subtract(triangle.get(0));
    Vertex P1Intersect = intersect.subtract(triangle.get(1));
    Vertex P2Intersect = intersect.subtract(triangle.get(2));
    
    if(edge1.crossProduct(P0Intersect).dotProduct(normal) < 0)  {
      return -1;
    }
    if(edge2.crossProduct(P1Intersect).dotProduct(normal) < 0)  {
      return -1;
    }
    if(edge3.crossProduct(P2Intersect).dotProduct(normal) < 0)  {
      return -1;
    }
    return tVal;
  }
  
  Vertex getIntersect(Ray r, float t)  {
    float xIntersect = r.getOrigin().v[0] + t * r.getVector().v[0];
    float yIntersect = r.getOrigin().v[1] + t * r.getVector().v[1];
    float zIntersect = r.getOrigin().v[2] + t * r.getVector().v[2];
    Vertex intersect = new Vertex(xIntersect, yIntersect, zIntersect);
    return intersect;
  }
  
  Vertex getIntersectNorm(Ray r, float t)  {
    Vertex edge1 = triangle.get(1).subtract(triangle.get(0));
    Vertex edge2 = triangle.get(2).subtract(triangle.get(1));
    Vertex normal = edge2.crossProduct(edge1);
    
    float normMag = sqrt(pow(normal.v[0], 2) + pow(normal.v[1], 2) + pow(normal.v[2], 2));
    normal.v[0] = normal.v[0] / normMag;
    normal.v[1] = normal.v[1] / normMag;
    normal.v[2] = normal.v[2] / normMag;
    
    return normal;
  }
}

class Sphere extends Shape  {
  float[] sphereDetails = new float[4];
  
  Sphere(float radius, float x, float y, float z)  {
    sphereDetails[0] = radius;
    sphereDetails[1] = x;
    sphereDetails[2] = y;
    sphereDetails[3] = z;
  } 
  
  float intersects(Ray r)  {
    float quadPlus, quadMinus;
    float a = r.getVector().dotProduct(r.getVector());
    Vertex spherePos = new Vertex(sphereDetails[1], sphereDetails[2], sphereDetails[3]);
    Vertex QO = new Vertex(r.getOrigin().subtract(spherePos).v[0], r.getOrigin().subtract(spherePos).v[1], r.getOrigin().subtract(spherePos).v[2]);
    float b = 2 * (QO.dotProduct(r.getVector()));
    float c = (QO.dotProduct(QO)) - pow(sphereDetails[0],2);
    float discriminant = b * b - 4 * a * c;
    
    if(discriminant < 0)  {
      return -1;
    }
    else  {
      quadPlus = (-1 * b + sqrt(discriminant)) / (2.0f * a);
      quadMinus = (-1 * b - sqrt(discriminant)) / (2.0f * a);
    }
    
    if(quadPlus > quadMinus && quadMinus > 0)  {
      return quadMinus;
    }
    return quadPlus;
  }
  
  Vertex getIntersect(Ray r, float t)  {
    float xIntersect = r.getOriginX() + t * r.getVectorX();
    float yIntersect = r.getOriginY() + t * r.getVectorY();
    float zIntersect = r.getOriginZ() + t * r.getVectorZ();
    Vertex P = new Vertex(xIntersect, yIntersect, zIntersect);
    
    return P;
  }
    
  Vertex getIntersectNorm(Ray r, float t)  {
    float xIntersect = r.getOriginX() + t * r.getVectorX();
    float yIntersect = r.getOriginY() + t * r.getVectorY();
    float zIntersect = r.getOriginZ() + t * r.getVectorZ();
    
    Vertex P = new Vertex(xIntersect, yIntersect, zIntersect);
    Vertex spherePos = new Vertex(sphereDetails[1], sphereDetails[2], sphereDetails[3]);
    Vertex sphereNorm = P.subtract(spherePos);
    
    float normMag = sqrt(pow(sphereNorm.v[0], 2) + pow(sphereNorm.v[1], 2) + pow(sphereNorm.v[2], 2));
    sphereNorm.v[0] = sphereNorm.v[0]/normMag;
    sphereNorm.v[1] = sphereNorm.v[1]/normMag;
    sphereNorm.v[2] = sphereNorm.v[2]/normMag;
    return sphereNorm;
  }
}

class Surface{
  float[] surfaceDetails;
  Surface(float cDr, float cDg, float cDb, float cAr, float cAg, float cAb, float cSr, float cSg, float cSb, float P, float kRefl){
    surfaceDetails = new float[11];
    surfaceDetails[0] = cDr;
    surfaceDetails[1] = cDg;
    surfaceDetails[2] = cDb;
    surfaceDetails[3] = cAr;
    surfaceDetails[4] = cAg;
    surfaceDetails[5] = cAb;
    surfaceDetails[6] = cSr;
    surfaceDetails[7] = cSg;
    surfaceDetails[8] = cSb;
    surfaceDetails[9] = P;
    surfaceDetails[10] = kRefl;  
  }
}

class Light{
  float[] light;
  Light(float x, float y, float z, float r, float g, float b){
    light = new float[6];
    light[0] = x;
    light[1] = y;
    light[2] = z;
    light[3] = r;
    light[4] = g;
    light[5] = b;
  }
}
    

///////////////////////////////////////////////////////////////////////
//
// Press key 1 to 9 and 0 to run different test cases.
//
///////////////////////////////////////////////////////////////////////
void keyPressed() {
  switch(key) {
    case '1':  gCurrentFile = new String("t0.cli"); interpreter(); break;
    case '2':  gCurrentFile = new String("t1.cli"); interpreter(); break;
    case '3':  gCurrentFile = new String("t2.cli"); interpreter(); break;
    case '4':  gCurrentFile = new String("t3.cli"); interpreter(); break;
    case '5':  gCurrentFile = new String("c0.cli"); interpreter(); break;
    case '6':  gCurrentFile = new String("c1.cli"); interpreter(); break;
    case '7':  gCurrentFile = new String("c2.cli"); interpreter(); break;
    case '8':  gCurrentFile = new String("c3.cli"); interpreter(); break;
    case '9':  gCurrentFile = new String("c4.cli"); interpreter(); break;
    case '0':  gCurrentFile = new String("c5.cli"); interpreter(); break;
  }
}

///////////////////////////////////////////////////////////////////////
//
//  Parser core. It parses the CLI file and processes it based on each 
//  token. Only "color", "rect", and "write" tokens are implemented. 
//  You should start from here and add more functionalities for your
//  ray tracer.
//
//  Note: Function "splitToken()" is only available in processing 1.25 
//       or higher.
//
///////////////////////////////////////////////////////////////////////
void interpreter() {
  bgR = 0;
  bgG = 0;
  bgB = 0;
  lights = new ArrayList<Light>();
  shapes = new ArrayList<Shape>();
  
  String str[] = loadStrings(gCurrentFile);
  if (str == null) println("Error! Failed to read the file.");
  for (int i=0; i<str.length; i++) {
    
    String[] token = splitTokens(str[i], " "); // Get a line and parse tokens.
    if (token.length == 0) continue; // Skip blank line.
    
    if (token[0].equals("fov")) {
      angle = radians(float(token[1]));
      k = tan(angle / 2);
    }
    else if (token[0].equals("background")) {
      bgR = float(token[1]);
      bgG = float(token[2]);
      bgB = float(token[3]);
    }
    else if (token[0].equals("light")) {
      float x = float(token[1]);
      float y = float(token[2]);
      float z = float(token[3]);
      float r = float(token[4]);
      float g = float(token[5]);
      float b = float(token[6]);
      Light light = new Light(x, y, z, r, g, b);
      lights.add(light);
      
    }
    else if (token[0].equals("surface")) {
      float cDr = float(token[1]);
      float cDg = float(token[2]);
      float cDb = float(token[3]);
      float cAr = float(token[4]);
      float cAg = float(token[5]);
      float cAb = float(token[6]);
      float cSr = float(token[7]);
      float cSg = float(token[8]);
      float cSb = float(token[9]);
      float P = float(token[10]);
      float kRefl = float(token[11]);
      surface = new Surface(cDr, cDg, cDb, cAr, cAg, cAb, cSr, cSg, cSb, P, kRefl);
    }    
    else if (token[0].equals("sphere")) {
      float r = float(token[1]);
      float x = float(token[2]);
      float y = float(token[3]);
      float z = float(token[4]);
      Sphere sphere = new Sphere(r, x, y, z);
      sphere.surface(surface);
      shapes.add(sphere);
    }
    else if (token[0].equals("begin")) {
      triangle = new Triangle();
    }
    else if (token[0].equals("vertex")) {
      float x = float(token[1]);
      float y = float(token[2]);
      float z = float(token[3]);
      triangle.vertex(x, y, z);
    }
    else if (token[0].equals("end")) {
      triangle.surface(surface);
      shapes.add(triangle);
    }
    else if (token[0].equals("color")) {
      float r = float(token[1]);
      float g = float(token[2]);
      float b = float(token[3]);
      fill(r, g, b);
    }
    else if (token[0].equals("rect")) {
      float x0 = float(token[1]);
      float y0 = float(token[2]);
      float x1 = float(token[3]);
      float y1 = float(token[4]);
      rect(x0, height-y1, x1-x0, y1-y0);
    }
    else if (token[0].equals("write")) {
      // you should render the scene here
      rayTracer();
      save(token[1]);  
    }
  }
}


///////////////////////////////////////////////////////////////////////
//
// Some initializations for the scene.
//
///////////////////////////////////////////////////////////////////////
void setup() {
  size(300, 300);  
  noStroke();
  colorMode(RGB, 1.0);
  background(0, 0, 0);
  interpreter();
}

///////////////////////////////////////////////////////////////////////
//
// Draw frames.  Should leave this empty.
//
///////////////////////////////////////////////////////////////////////
void draw() {
}

void rayTracer()  {
  color c;
  Ray tempRay = new Ray();
  loadPixels();
  
  for(int i = 0; i < width; i++)  {
    for(int j = 0; j < height; j++)  {
      float xPrime = (i - width / 2) * (2 * k / width);
      float yPrime = -(j - height / 2) * (2 * k / height);
      float vMag = sqrt((pow(xPrime, 2) + pow(yPrime, 2) + 1));
      tempRay.origin(0, 0, 0);
      tempRay.vector(xPrime/vMag, yPrime/vMag, -1/vMag);
      
      float minTime = MAX_FLOAT;
      boolean hit = false;
      Shape oColor = null;
      
      for(Shape s : shapes)  {
        float time = s.intersects(tempRay);
        if(time > 0 && time < minTime)  {
          minTime = time;  
          oColor = s;
          hit = true;
        }
      }
    
      if(hit)  {
        Vertex baseColor = new Vertex(oColor.shapeDetails.surfaceDetails[3], oColor.shapeDetails.surfaceDetails[4], oColor.shapeDetails.surfaceDetails[5]);
        for(Light l : lights)  {
          Vertex intersectPoint = oColor.getIntersect(tempRay,minTime);
          Vertex lightPos = new Vertex(l.light[0], l.light[1], l.light[2]);
          float lightDist = intersectPoint.euclideanDistance(lightPos);
          Vertex surfaceColor = new Vertex(oColor.shapeDetails.surfaceDetails[0], oColor.shapeDetails.surfaceDetails[1], oColor.shapeDetails.surfaceDetails[2]);
          
          Vertex L = lightPos.subtract(oColor.getIntersect(tempRay, minTime));
          float Lmag = sqrt(pow(L.v[0], 2) + pow(L.v[1], 2) + pow(L.v[2], 2));
          L.v[0] = L.v[0] / Lmag;
          L.v[1] = L.v[1] / Lmag;
          L.v[2] = L.v[2] / Lmag;
          
          Ray lRay = new Ray();
          lRay.origin = new Vertex(intersectPoint.v[0], intersectPoint.v[1], intersectPoint.v[2]);
          lRay.vector = new Vertex(L.v[0], L.v[1], L.v[2]);
          float shade = MAX_FLOAT;
          Shape sIntersect = null;
          for(Shape sh : shapes)  {
            float curr = sh.intersects(lRay);
            if(curr > 0 && curr < shade && sh != oColor)  {
              sIntersect = sh;
              shade = curr;
            }
          }
          float lightInt = 0;
          if(sIntersect != null)  {
            Vertex intersect2 = sIntersect.getIntersect(lRay,shade);
            lightInt = intersect2.euclideanDistance(lightPos);
          }
          
          Vertex origin = new Vertex(0,0,0);
          Vertex surfaceNorm = (oColor.getIntersectNorm(tempRay, minTime));
          float surfaceNormMag = sqrt(pow(surfaceNorm.v[0], 2) + pow(surfaceNorm.v[1], 2) + pow(surfaceNorm.v[2], 2));
          surfaceNorm.v[0] = surfaceNorm.v[0] / surfaceNormMag;
          surfaceNorm.v[1] = surfaceNorm.v[1] / surfaceNormMag;
          surfaceNorm.v[2] = surfaceNorm.v[2] / surfaceNormMag;
          
          Vertex V = origin.subtract(oColor.getIntersect(tempRay, minTime));
          float VMag = sqrt(pow(V.v[0], 2) + pow(V.v[1], 2) + pow(V.v[2], 2));
          V.v[0] = V.v[0] / VMag;
          V.v[1] = V.v[1] / VMag;
          V.v[2] = V.v[2] / VMag;

          float diffuse = max(0, surfaceNorm.dotProduct(L));
          float NL = 2.0 * (surfaceNorm.dotProduct(L));
          
          surfaceNorm.v[0] = surfaceNorm.v[0] * NL;
          surfaceNorm.v[1] = surfaceNorm.v[1] * NL;
          surfaceNorm.v[2] = surfaceNorm.v[2] * NL;
          
          Vertex R = surfaceNorm.subtract(L);
          float sp2 = max(0, V.dotProduct(R));
          float specular = pow(sp2, oColor.shapeDetails.surfaceDetails[9]);
          
          Vertex kS = new Vertex(oColor.shapeDetails.surfaceDetails[6], oColor.shapeDetails.surfaceDetails[7], oColor.shapeDetails.surfaceDetails[8]);
          kS.v[0] = kS.v[0] * specular;
          kS.v[1] = kS.v[1] * specular;
          kS.v[2] = kS.v[2] * specular;
          
          surfaceColor.v[0] = surfaceColor.v[0] * diffuse;
          surfaceColor.v[1] = surfaceColor.v[1] * diffuse;
          surfaceColor.v[2] = surfaceColor.v[2] * diffuse;
          
          Vertex addColor = surfaceColor.add(kS);
          Vertex lColor = new Vertex(l.light[3], l.light[4], l.light[5]);
          
          addColor.v[0] = addColor.v[0] * l.light[3];
          addColor.v[1] = addColor.v[1] * l.light[4];
          addColor.v[2] = addColor.v[2] * l.light[5];
          
          if(sIntersect == null)  {
            baseColor = baseColor.add(addColor);
          }
        }
        c = color(baseColor.v[0], baseColor.v[1], baseColor.v[2]);
      }
      else  {
        c = color(bgR, bgG, bgB);
      }
      pixels[j * width + i] = c; 
    }
  }
  updatePixels();
}
