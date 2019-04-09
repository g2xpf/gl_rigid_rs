use crate::circle;
use crate::rectangle;

#[derive(Debug)]
pub enum E {
    Circle(circle::Circle),
    Rectangle(rectangle::Rectangle),
}

impl From<circle::Circle> for E {
    fn from(c: circle::Circle) -> E {
        E::Circle(c)
    }
}

impl From<rectangle::Rectangle> for E {
    fn from(r: rectangle::Rectangle) -> E {
        E::Rectangle(r)
    }
}

pub struct GlobalUniform {
    pub dt: f32,
    pub gravity: [f32; 2],
    pub scale: f32,
    pub resolution: [f32; 2],
}

#[derive(Debug)]
pub struct Collision {
    pub lhs: usize,
    pub rhs: usize,
    pub normal: [f32; 2],
    pub depth: f32,
    pub collision_point: [f32; 2],
}
