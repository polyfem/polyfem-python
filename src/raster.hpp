#pragma once

#include <Eigen/Core>

#include <functional>
#include <vector>

namespace renderer
{
	struct Material
	{
		Eigen::Vector3d diffuse_color;
		Eigen::Vector3d specular_color;
		double specular_exponent;
	};
	class VertexAttributes
	{
	public:
		VertexAttributes(double x = 0, double y = 0, double z = 0, double w = 1)
		{
			position << x, y, z, w;
		}

		static VertexAttributes interpolate(
			const VertexAttributes &a,
			const VertexAttributes &b,
			const VertexAttributes &c,
			const double alpha,
			const double beta,
			const double gamma)
		{
			VertexAttributes r;
			r.position = alpha * a.position + beta * b.position + gamma * c.position;
			r.color = alpha * a.color + beta * b.color + gamma * c.color;
			return r;
		}

		Eigen::Vector4d position;
		Eigen::Vector3d normal;
		Eigen::Vector3d color;
		Material material;
	};

	class FragmentAttributes
	{
	public:
		FragmentAttributes(double r = 0, double g = 0, double b = 0, double a = 1)
		{
			color << r, g, b, a;
		}
		double depth;
		Eigen::Vector4d color;
	};

	class FrameBufferAttributes
	{
	public:
		FrameBufferAttributes(uint8_t r = 0, uint8_t g = 0, uint8_t b = 0, uint8_t a = 255)
		{
			color << r, g, b, a;
			depth = 2;
		}
		double depth;
		Eigen::Matrix<uint8_t, 4, 1> color;
	};

	class UniformAttributes
	{
	public:
		Eigen::Matrix4d M_cam;

		Eigen::Matrix4d M_orth;
		Eigen::Matrix4d M;
	};

	class Program
	{
	public:
		std::function<VertexAttributes(const VertexAttributes &, const UniformAttributes &)> VertexShader;
		std::function<FragmentAttributes(const VertexAttributes &, const UniformAttributes &)> FragmentShader;
		std::function<FrameBufferAttributes(const FragmentAttributes &, const FrameBufferAttributes &)> BlendingShader;
	};

	void rasterize_triangle(const Program &program, const UniformAttributes &uniform, const VertexAttributes &v1, const VertexAttributes &v2, const VertexAttributes &v3, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer);
	void rasterize_triangles(const Program &program, const UniformAttributes &uniform, const std::vector<VertexAttributes> &vertices, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer);

	void rasterize_line(const Program &program, const UniformAttributes &uniform, const VertexAttributes &v1, const VertexAttributes &v2, double line_thickness, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer);
	void rasterize_lines(const Program &program, const UniformAttributes &uniform, const std::vector<VertexAttributes> &vertices, double line_thickness, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer);

	void framebuffer_to_uint8(const Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer, std::vector<uint8_t> &image);
}
