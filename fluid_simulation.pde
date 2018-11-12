// Implementation of algorithm described in Real-Time Fluid Dynamics for Games by Jos Stam http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf

float visc = 0.0f;
float diff = 0.0f;

int N = 250;
int size = (N + 2) * (N + 2);
float[] u = new float[size];
float[] v = new float[size];
float[] u_prev = new float[size];
float[] v_prev = new float[size];
float[] dens = new float[size];
float[] dens_prev = new float[size];

float dt = 0.4f;
float densityAdd = 50.0f;
float forceX = 0.0f;
float forceY = -7.0f;
float circleRadius = 1.0f;

void setup()
{
  size(252, 252);
}

int IX(int i, int j)
{
  return i + (N + 2) * j; 
}

void addSource(float[] x, float s[], float dt)
{
  for(int index = 0; index < size; index++)
  {
    x[index] += dt * s[index];
  }
}

void diffuse(int N, int b, float[] x, float[] x0, float diff, float dt)
{
  int i, j, k;
  float a = dt * diff * N * N;
  for(k = 0; k < 20; k++)
  {
    for(i = 1; i <= N; i++)
    {
      for(j = 1; j <= N; j++)
      {
        x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / (1 + 4 * a);
      }
    }
    setBounds(N, b, x);
  }
}

void advect(int N, int b, float[] d, float[] d0, float[] u, float v[], float dt)
{
  int i, j, i0, j0, i1, j1;
  float x, y, s0, t0, s1, t1;
  float dt0 = dt * N;
  
  for(i = 1; i <= N; i++)
  {
    for(j = 1; j <= N; j++)
    {
      x = i - dt0 * u[IX(i, j)];
      y = j - dt0 * v[IX(i, j)];
      
      if (x < 0.5)
      {
        x = 0.5;
      }
      if (x > N + 0.5)
      {
        x = N + 0.5;
      }
      
      i0 = (int)x;
      i1 = i0 + 1;
      
      if (y < 0.5)
      {
        y = 0.5;
      }
      if (y > N + 0.5)
      {
        y = N + 0.5;
      }
      j0 = (int)y;
      j1 = j0 + 1;
      
      s1 = x-i0;
      s0 = 1-s1;
      t1 = y-j0;
      t0 = 1-t1;
      d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) + s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
    }
  }
  
  setBounds(N, b, d);
}

void densityStep(int N, float[] x, float[] x0, float[] u, float[] v, float diff, float dt)
{
  float[] tmp;
  addSource(x, x0, dt);
  tmp = x0; x0 = x; x = tmp; // Swap
  diffuse(N, 0, x, x0, diff, dt);
  tmp = x0; x0 = x; x = tmp; // Swap Back
  advect(N, 0, x, x0, u, v, dt);
}

void project(int N, float[] u, float[] v, float[] p, float[] div)
{
  int i, j, k;
  float h = 1.0f / N;
  
  for (i = 1; i <= N; i++)
  {
    for (j = 1; j <= N; j++)
    {
      div[IX(i, j)] = -0.5 * h * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]);
      p[IX(i, j)] = 0;
    }
  }
  
  setBounds(N, 0, div);
  setBounds(N, 0, p);
  
  for(k = 0; k < 20; k++)
  {
    for(i = 1 ;i <= N; i++)
    {
      for(j = 1; j <= N; j++)
      {
        p[IX(i, j)] = (div[IX(i, j)] + p[IX(i - 1, j)] + p[IX(i + 1, j)] + p[IX(i, j - 1)] + p[IX(i, j + 1)]) / 4.0f;
      }
    }
    setBounds(N, 0, p);
  }
  
  for(i = 1; i <= N; i++)
  {
    for(j = 1; j <= N; j++)
    {
      u[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
      v[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
    }
  }
  
  setBounds(N, 1, u);
  setBounds(N, 2, v);
}

void velocityStep(int N, float[] u, float[] v, float[] u0, float[] v0, float visc, float dt)
{
  float[] tmp;
  
  addSource(u, u0, dt);
  addSource(v, v0, dt);
  tmp = u0; u0 = u; u = tmp; // Swap
  diffuse(N, 1, u, u0, visc, dt);
  tmp = v0; v0 = v; v = tmp; // Swap
  diffuse(N, 2, v, v0, visc, dt);
  project(N, u, v, u0, v0);
  tmp = u0; u0 = u; u = tmp; // Swap Back
  tmp = v0; v0 = v; v = tmp; // Swap Back
  advect(N, 1, u, u0, u0, v0, dt);
  advect(N, 2, v, v0, u0, v0, dt);
  project(N, u, v, u0, v0);
}

void setBounds (int N, int b, float[] x)
{
  for(int i = 1; i <= N; i++)
  {
    x[IX(0, i)] = b==1 ? -x[IX(1, i)] : x[IX(1, i)];
    x[IX(N + 1, i)] = b==1 ? -x[IX(N, i)] : x[IX(N, i)];
    x[IX(i, 0)] = b==2 ? -x[IX(i, 1)] : x[IX(i, 1)];
    x[IX(i, N + 1)] = b==2 ? -x[IX(i, N)] : x[IX(i, N)];
  }
  
  x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
  x[IX(0, N + 1)] = 0.5 * (x[IX(1,N + 1)] + x[IX(0, N)]);
  x[IX(N + 1, 0)] = 0.5 * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
  x[IX(N + 1, N + 1)] = 0.5 * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

void updateSimulation(float dt)
{
  velocityStep(N, u, v, u_prev, v_prev, visc, dt);
  densityStep(N, dens, dens_prev, u, v, diff, dt);
}

void drawDensity()
{
  background(0);
  loadPixels();
  
  float density;
  for(int index = 0; index < size; index++)
  {
    density = dens[index] * 255.0f;
    pixels[index] = color(density, density, density);
  }
  updatePixels();
}

void drawCircle(float dt, float val, float radius)
{
  // create B&W mask and apply density where pixels > 0.0f
  background(0);
  noStroke();
  ellipse(mouseX, mouseY, radius, radius);
  loadPixels();
  
  float pixelVal;
  for(int index = 0; index < size; index++)
  {
    pixelVal = red(pixels[index]);
    
    if(pixelVal > 0.0f)
    {
      dens_prev[index] = val;
      u_prev[index] += dt * forceX;
      v_prev[index] += dt * forceY;
    }
  }
}

void draw()
{  
  for(int index = 0; index < size; index++)
  {
    dens_prev[index] = 0;
    u_prev[index] = 0;
    v_prev[index] = 0;
  }
  
  if (mousePressed == true)
  {
    if(mouseButton == LEFT)
    {
      drawCircle(dt, densityAdd, circleRadius);
    }
  }
  
  updateSimulation(dt);
  drawDensity();
}
