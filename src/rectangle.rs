use glium::{index, vertex};

static RECTANGLE_VS_SOURCE: &'static str = include_str!("rectangle.vert");
static RECTANGLE_FS_SOURCE: &'static str = include_str!("rectangle.frag");

#[derive(Copy, Clone)]
pub struct Vertex {
    pub coord: [f32; 2],
}
implement_vertex!(Vertex, coord);

#[derive(Debug)]
pub struct Rectangle {
    pub r: [f32; 2],
    pub v: [f32; 2],
    pub a: [f32; 2],
    pub wh: [f32; 2],
    pub omega: f32,
    pub angle: f32,
    pub is_static: bool,
    pub mass: f32,
    pub inv_mass: f32,
    pub inertia: f32,
    pub inv_inertia: f32,

    pub visible: bool,
}

pub struct RectangleBrander<'a> {
    pub program: glium::Program,
    pub index_buffer: index::IndexBuffer<u32>,
    pub vertex_buffer: vertex::VertexBuffer<Vertex>,
    pub draw_parameter: glium::DrawParameters<'a>,
}

impl Rectangle {
    pub fn new(r: [f32; 2], v: [f32; 2], wh: [f32; 2], is_static: bool) -> Self {
        let mass = if is_static { wh[0] * wh[1] } else { 0.0 };
        let inertia = 1.0 / 12.0 * mass * (wh[0] * wh[0] + wh[1] * wh[1]);
        Rectangle {
            r: r,
            v: v,
            a: [0.0, 0.0],
            wh: wh,
            omega: 0.0,
            angle: 0.0,
            is_static: is_static,
            mass: mass,
            inv_mass: if is_static { 0.0 } else { 1.0 / mass },
            inertia: inertia,
            inv_inertia: if is_static { 0.0 } else { 1.0 / inertia },

            visible: true,
        }
    }
}

impl<'a> RectangleBrander<'a> {
    pub fn new(display: &glium::Display, draw_parameter: glium::DrawParameters<'a>) -> Self {
        RectangleBrander {
            program: glium::Program::from_source(
                display,
                RECTANGLE_VS_SOURCE,
                RECTANGLE_FS_SOURCE,
                None,
            )
            .unwrap(),
            vertex_buffer: glium::VertexBuffer::new(
                display,
                &[
                    Vertex { coord: [1.0, 1.0] },
                    Vertex { coord: [-1.0, 1.0] },
                    Vertex { coord: [1.0, -1.0] },
                    Vertex {
                        coord: [-1.0, -1.0],
                    },
                ],
            )
            .unwrap(),
            index_buffer: glium::IndexBuffer::new(
                display,
                index::PrimitiveType::TrianglesList,
                &[0, 1, 2, 1, 2, 3],
            )
            .unwrap(),
            draw_parameter: draw_parameter,
        }
    }
}
