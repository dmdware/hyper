
#include "../sys/includes.h"
#include "../sys/texture.h"
#include "../render/ms3d.h"


void rotate(float* v, float rad, float x, float y, float z, float* newv)
{
	float cosTheta = (float)cos(rad);
	float sinTheta = (float)sin(rad);

	newv[0] = (cosTheta + (1 - cosTheta) * x * x)		* v[0];
	newv[0] += ((1 - cosTheta) * x * y - z * sinTheta)	* v[1];
	newv[0] += ((1 - cosTheta) * x * z + y * sinTheta)	* v[2];

	newv[1] = ((1 - cosTheta) * x * y + z * sinTheta)	* v[0];
	newv[1] += (cosTheta + (1 - cosTheta) * y * y)		* v[1];
	newv[1] += ((1 - cosTheta) * y * z - x * sinTheta)	* v[2];

	newv[2] = ((1 - cosTheta) * x * z - y * sinTheta)	* v[0];
	newv[2] += ((1 - cosTheta) * y * z + x * sinTheta)	* v[1];
	newv[2] += (cosTheta + (1 - cosTheta) * z * z)		* v[2];
}

float dot(Vec3f vVector1, Vec3f vVector2)
{
	return ((vVector1.x * vVector2.x) + (vVector1.y * vVector2.y) + (vVector1.z * vVector2.z));
}

void makeplane(float* d, float* inp, float* n)
{
	float mag = sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	n[0] /= mag;
	n[1] /= mag;
	n[2] /= mag;
	*d = -dot(*(Vec3f*)n, *(Vec3f*)inp);
}

bool pointonorbehindplane(Vec3f point, Vec3f planenormal, float planed, float epsilon)
{
	float result = point.x*planenormal.x + point.y*planenormal.y + point.z*planenormal.z + planed;

	if (result <= epsilon)
		return true;

	return false;
}

Vec3f pop(Vec3f n, float d)
{
	Vec3f p;
	p.x = n.x * -d;
	p.y = n.y * -d;
	p.z = n.z * -d;
	//printf("pop %f,%f,%f\r\n", p.x, p.y, p.z);
	return p;
}

Vec3f tov(Vec3f a, Vec3f b)
{
	Vec3f c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
	return c;
}

float planedistance(Vec3f normal, float dist, Vec3f point)
{
#if 0
	Vec3f ppoint = pop(normal, dist);
	float pd = dot(normal, (tov(point, ppoint)));
#else
	Vec3f pp = pop(normal, dist);
	Vec3f v2 = tov(point, pp);
	float pd = -dot(v2, normal);
#endif
	printf("pd %f\r\n", pd);
	return pd;
}

float sign(Vec3f p1, Vec3f p2, Vec3f p3)
{
	return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

bool pointintriangle(Vec3f pt, Vec3f v1, Vec3f v2, Vec3f v3)
{
	bool b1, b2, b3;

	b1 = sign(pt, v1, v2) < 0.0f;
	b2 = sign(pt, v2, v3) < 0.0f;
	b3 = sign(pt, v3, v1) < 0.0f;

	return ((b1 == b2) && (b2 == b3));
}

float mag(Vec3f v)
{
	return (float)sqrtf((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
}

Vec3f norm(Vec3f n)
{
	float m;

	m = mag(n);

	n.x /= m;
	n.y /= m;
	n.z /= m;

	return n;
}

Vec3f cross(Vec3f v1, Vec3f v2)
{
	Vec3f n;

	n.x = ((v1.y * v2.z) - (v1.z * v2.y));
	n.y = ((v1.z * v2.x) - (v1.x * v2.z));
	n.z = ((v1.x * v2.y) - (v1.y * v2.x));

	return n;
}

Vec3f tnorm(Vec3f* t)
{
	Vec3f v1;
	Vec3f v2;
	Vec3f n;
	v1 = tov(t[2], t[0]);
	v2 = tov(t[1], t[0]);
	n = cross(v1, v2);
	n = norm(n);
	return n;
}

Vec3f projvecontopl(Vec3f v, Vec3f pn, Vec3f pp)
{
#if 0
	v.x -= pp.x;
	v.y -= pp.y;
	v.z -= pp.z;
	float d = dot(v, pn);
	Vec3f r;
	//assume pn normalized
	r.x = (v.x - pn.x * d);
	r.y = (v.y - pn.y * d);
	r.z = (v.z - pn.z * d);
	v.x += pp.x;
	v.y += pp.y;
	v.z += pp.z;
#else
	Vec3f v2 = tov(v, pp);
	float d = dot(v2, pn);
	Vec3f r;
	r.x = (v.x - pn.x * d);
	r.y = (v.y - pn.y * d);
	r.z = (v.z - pn.z * d);
#endif
	return r;
}

void main()
{
	MS3DModel m;
	unsigned int test;
	char c[3];
	m.load("asd.ms3d", test, test, test, test, true);

#define WX	256
#define WY	512

	unsigned char rgb[WX*WY * 3];
	float depth[WX*WY];

	for (int i = 0; i < WX*WY; ++i)
	{
		rgb[i * 3 + 0] = 255;
		rgb[i * 3 + 1] = 255;
		rgb[i * 3 + 2] = 255;
		depth[i] = 1;
	}

#define NEAR		1
#define FAR			110
#define ASPECT		((float)WX/(float)WY)
#define FOV			45
#define POSX		-100
#define POSY		150
#define POSZ		-150
#define VIEWX		0
#define VIEWY		0
#define VIEWZ		0
#define UPX			0
#define UPY			1
#define UPZ			0
#define RIGHTX		-150
#define RIGHTY		0
#define RIGHTZ		100
#define DISTSCALE	1
#define OFFX		5
#define OFFY		0
#define OFFZ		0

	Vec3f up;
	Vec3f right;
	Vec3f viewdir;
	Vec3f cpos;
	Vec3f viewpos;
	cpos.x = POSX;
	cpos.y = POSY;
	cpos.z = POSZ;
	viewpos.x = VIEWX;
	viewpos.y = VIEWY;
	viewpos.z = VIEWZ;
	viewdir = tov(viewpos, cpos);
	up.x = UPX;
	up.y = UPY;
	up.z = UPZ;
	right.x = RIGHTX;
	right.y = RIGHTY;
	right.z = RIGHTZ;
	right = norm(right);
	up = norm(cross(right, viewdir));

	for (int i = 0; i < m.m_numVertices; ++i)
	{
		float worldv[3];
		worldv[0] = m.m_pVertices[i].m_location[0] + OFFX;
		worldv[1] = m.m_pVertices[i].m_location[1] + OFFY;
		worldv[2] = m.m_pVertices[i].m_location[2] + OFFZ;

		float clipv[3];
		float offv[3];
		offv[0] = worldv[0] - cpos.x;
		offv[1] = worldv[1] - cpos.y;
		offv[2] = worldv[2] - cpos.z;
		float depthf = sqrtf(offv[0] * offv[0] + offv[1] * offv[1] + offv[2] * offv[2]);

#if 01
		float sideplaned;
		Vec3f sideplanen;
		float vertplaned;
		Vec3f vertplanen;
		sideplanen.x = right.x;
		sideplanen.y = right.y;
		sideplanen.z = right.z;
		sideplanen = norm(sideplanen);
		vertplanen.x = up.x;
		vertplanen.y = up.y;
		vertplanen.z = up.z;
		vertplanen = norm(vertplanen);
		float pos[3];
		pos[0] = cpos.x;
		pos[1] = cpos.y;
		pos[2] = cpos.z;
		makeplane(&sideplaned, pos, (float*)&sideplanen);
		makeplane(&vertplaned, pos, (float*)&vertplanen);
#if 0
		float sideouter = (depthf / DISTSCALE) / 1.0;
		float vertouter = (depthf / DISTSCALE) / 1.0f;
		//sideouter *= sideouter;
		//vertouter *= vertouter;
		sideouter = (pow(0.0f + sideouter, 2.0f) - 0.0f) * ASPECT / 4000000.0f;
		vertouter = (pow(0.0f + vertouter, 2.0f) - 0.0f) / 4000000.0f;
#elif 0
		float sideouter = ((depthf) * ASPECT / 5.0f + pow(depthf / 0.01f, 5.0f) * ASPECT) / 200000000000000000000.0f;
		float vertouter = ((depthf) / 5.0f + pow(depthf / 0.01f, 5.0f)) / 200000000000000000000.0f;
#elif 01
		float sideouter = ((depthf)* ASPECT / 5.0f + pow(depthf / 0.01f, 6.0f) * ASPECT) / 5000000000000000000000000.0f;
		float vertouter = ((depthf) / 5.0f + pow(depthf / 0.01f, 6.0f)) / 5000000000000000000000000.0f;
#elif 0
		float sideouter = (depthf)* ASPECT / 5.0f;
		float vertouter = (depthf) / 5.0f;
#endif

		float sideways;
		float vertways;

		if (pointonorbehindplane(*(Vec3f*)worldv, sideplanen, sideplaned, 0))
		{
			sideways = planedistance(sideplanen, sideplaned, *(Vec3f*)worldv);
		}
		else
		{
			sideways = planedistance(sideplanen, sideplaned, *(Vec3f*)worldv);
		}

		if (pointonorbehindplane(*(Vec3f*)worldv, vertplanen, vertplaned, 0))
		{
			vertways = planedistance(vertplanen, vertplaned, *(Vec3f*)worldv);
		}
		else
		{
			vertways = planedistance(vertplanen, vertplaned, *(Vec3f*)worldv);
		}

		depthf /= FAR;

		clipv[0] = ((sideways / sideouter) + 1.0f) / 2.0f;
		clipv[1] = ( (vertways / vertouter) + 1.0f ) / 2.0f;
		clipv[2] = depthf;
#else
		//clipv[0] = 0.5f;
		//clipv[1] = 0.5f;
		//depthf = 0.5f;

		clipv[0] = (float)(rand() % 1000) / 1000.0f;
		clipv[1] = (float)(rand() % 1000) / 1000.0f;
		clipv[2] = (float)(rand() % 1000) / 1000.0f;

#endif

		m.m_pVertices[i].m_location[0] = clipv[0];
		m.m_pVertices[i].m_location[1] = clipv[1];
		m.m_pVertices[i].m_location[2] = clipv[2];
	}


	for (int i = 0; i < m.m_numTriangles; ++i)
	{
		MS3DModel::Triangle* t = &m.m_pTriangles[i];
		float vmin[3];
		float vmax[3];

		for (int j = 0; j < 3; ++j)
		{
			float *v = m.m_pVertices[t->m_vertexIndices[j]].m_location;

			if (v[0] < vmin[0] || j == 0)
				vmin[0] = v[0];
			if (v[1] < vmin[1] || j == 0)
				vmin[1] = v[1];
			if (v[2] < vmin[2] || j == 0)
				vmin[2] = v[2];
			if (v[0] > vmax[0] || j == 0)
				vmax[0] = v[0];
			if (v[1] > vmax[1] || j == 0)
				vmax[1] = v[1];
			if (v[2] > vmax[2] || j == 0)
				vmax[2] = v[2];
		}

		float tplaned;
		Vec3f tplanen;
		
		if (vmin[0] < 0)
			vmin[0] = 0;
		if (vmin[1] < 0)
			vmin[1] = 0;
		if (vmax[0] > 1)
			vmax[0] = 1;
		if (vmax[1] > 1)
			vmax[1] = 1;

		makeplane(&tplaned, m.m_pVertices[t->m_vertexIndices[0]].m_location, (float*)&tplanen);

		for (float sx = vmin[0]; sx <= vmax[0]; sx += 1.0f / WX)
		{
			for (float sy = vmin[1]; sy <= vmax[1]; sy += 1.0f / WY)
			{
				Vec3f sp;
				sp.x = sx;
				sp.y = sy;

				if (sp.x < 0 || sp.x >= 1 || sp.y < 0 || sp.y >= 1)
					continue;

				if (!pointintriangle(sp,
					*(Vec3f*)m.m_pVertices[t->m_vertexIndices[0]].m_location,
					*(Vec3f*)m.m_pVertices[t->m_vertexIndices[1]].m_location,
					*(Vec3f*)m.m_pVertices[t->m_vertexIndices[2]].m_location))
					continue;

				//printf("pit %f,%f,%f\r\n",
				//	m.m_pVertices[t->m_vertexIndices[0]].m_location[0],
				//	m.m_pVertices[t->m_vertexIndices[0]].m_location[1],
				//	m.m_pVertices[t->m_vertexIndices[0]].m_location[2]);

				sp.z = 0;
				Vec3f pp = pop(tplanen, tplaned);
				Vec3f dp = projvecontopl(sp, tplanen, pp);
				sp.z = dp.z;

				//printf("dp1 %f\r\n", sp.z);

				int pixi = ((int)(WX*sp.x) + WX*(int)(WY*sp.y));

				if (sp.z > depth[pixi] || sp.z < 0 || sp.z >= 1)
					continue;

				//printf("dp2\r\n");

				depth[pixi] = sp.z;

				rgb[pixi * 3 + 0] = 255 * sp.z;
				rgb[pixi * 3 + 1] = 255 * sp.z;
				rgb[pixi * 3 + 2] = 255 * sp.z;

				//rgb[pixi * 3 + 0] = 0;
				//rgb[pixi * 3 + 1] = 0;
				//rgb[pixi * 3 + 2] = 0;
			}
		}
	}

	//rgb[3 * (WX / 2 + WX * WY / 2) + 0] = 0;
	//rgb[3 * (WX / 2 + WX * WY / 2) + 1] = 0;
	//rgb[3 * (WX / 2 + WX * WY / 2) + 2] = 0;

	savepng("out.png", (unsigned char*)rgb, WX, WY, 3);

	printf("done\r\n");
	fgets(c, 2, stdin);
}







