
#include "AreaLight.hpp"


namespace FW {


void AreaLight::draw(const Mat4f& worldToCamera, const Mat4f& projection) {
    glUseProgram(0);
    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf((float*)&projection);
    glMatrixMode(GL_MODELVIEW);
    Mat4f S = Mat4f::scale(Vec3f(m_size,1));
    Mat4f M = worldToCamera *m_xform * S;
    glLoadMatrixf((float*)&M);
    glBegin(GL_TRIANGLES);
    glColor3fv( &m_E.x );
    glVertex3f(1,1,0); glVertex3f(1,-1,0); glVertex3f( -1,-1,0 );
    glVertex3f(1,1,0); glVertex3f( -1,-1,0 ); glVertex3f(-1,1,0); 
    glEnd();
}

void AreaLight::sample(float& pdf, Vec3f& p, int baseX, int baseY, int samplingIndex, Random& rnd) {
    // You should draw a random point on the light source and evaluate the PDF.
    // Store the results in "pdf" and "p".
    // 
    // Note: The "size" member is _one half_ the diagonal of the light source.
    // That is, when you map the square [-1,1]^2 through the scaling matrix
    // 
    // S = ( size.x    0   )
    //     (   0    size.y )
    // 
    // the result is the desired light source quad (see draw() function above).
    // This means the total area of the light source is 4*size.x*size.y.
    // This has implications for the computation of the PDF.

    // For extra credit, implement QMC sampling using some suitable sequence.
    // Use the "base" input for controlling the progression of the sequence from
    // the outside. If you only implement purely random sampling, "base" is not required.

    // Inefficient QMC based on Halton sequence
    // Using different base for x and y to avoid correlation
    float x = 0;
    float y = 0;
    float fX = 1.0 / baseX;
    float fY = 1.0 / baseY;
    int i = samplingIndex;

    while (i > 0) {
		x += fX * (i % baseX);
        i /= baseX;
		fX /= baseX;
	}
    i = samplingIndex;
    while (i > 0) {
        y += fY * (i % baseY);
		i /= baseY;
        fY /= baseY;
    }

    Mat4f objectToWorld = m_xform * Mat4f::scale(Vec3f(m_size, 1));
    Vec4f rndPointInObjectSpace = Vec4f(2 * x - 1, 2 * y - 1, 0, 1);
    p = (objectToWorld * rndPointInObjectSpace).getXYZ();
    pdf = 1.0f / (4 * m_size.x * m_size.y);
}


} // namespace FW
