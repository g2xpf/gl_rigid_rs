use glium::{index, vertex};

static CIRCLE_VS_SOURCE: &'static str = include_str!("circle.vert");
static CIRCLE_FS_SOURCE: &'static str = include_str!("circle.frag");

#[derive(Copy, Clone)]
pub struct Vertex {
    pub coord: [f32; 2],
}
implement_vertex!(Vertex, coord);

#[derive(Debug)]
pub struct Circle {
    pub r: [f32; 2],
    pub v: [f32; 2],
    pub a: [f32; 2],
    pub is_static: bool,
    pub omega: f32,
    pub angle: f32,
    pub radius: f32,
    pub mass: f32,
    pub inv_mass: f32,
    pub inertia: f32,
    pub inv_inertia: f32,

    pub visible: bool,
}

pub struct CircleBrander<'a> {
    pub program: glium::Program,
    pub index_buffer: index::IndexBuffer<u32>,
    pub vertex_buffer: vertex::VertexBuffer<Vertex>,
    pub draw_parameter: glium::DrawParameters<'a>,
}

impl Circle {
    pub fn new(radius: f32, r: [f32; 2], v: [f32; 2], is_static: bool) -> Self {
        let mass = if is_static {
            0.0
        } else {
            1.0 * radius * radius * std::f32::consts::PI
        };
        let inertia = 0.5 * mass * radius * radius;
        Circle {
            r: r,
            v: v,
            a: [0.0, 0.0],
            omega: 0.0,
            angle: 0.0,
            radius: radius,
            is_static: is_static,
            mass: mass,
            inv_mass: if is_static { 0.0 } else { 1.0 / mass },
            inertia: inertia,
            inv_inertia: if is_static { 0.0 } else { 1.0 / inertia },

            visible: true,
        }
    }
}

impl<'a> CircleBrander<'a> {
    pub fn new(display: &glium::Display, draw_parameter: glium::DrawParameters<'a>) -> Self {
        let split = 100;
        CircleBrander {
            program: glium::Program::from_source(display, CIRCLE_VS_SOURCE, CIRCLE_FS_SOURCE, None)
                .unwrap(),
            vertex_buffer: {
                let mut v = Vec::new();
                for i in 0..split {
                    let theta = (i as f32) * 2.0 * std::f32::consts::PI / (split as f32);
                    v.push(Vertex {
                        coord: [theta.cos(), theta.sin()],
                    })
                }
                v.push(Vertex { coord: [0.0, 0.0] });
                glium::VertexBuffer::new(display, &v[..]).unwrap()
            },
            index_buffer: {
                let mut v = Vec::new();
                for i in 0..split {
                    v.push(i);
                    v.push((i + 1) % split);
                    v.push(split);
                }
                glium::IndexBuffer::new(display, index::PrimitiveType::TrianglesList, &v[..])
                    .unwrap()
            },
            draw_parameter: draw_parameter,
        }
    }
}
