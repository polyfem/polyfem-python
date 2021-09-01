#include "raster.hpp"



#include <iostream>

namespace renderer
{
	void rasterize_triangle(const Program &program, const UniformAttributes &uniform, const VertexAttributes &v1, const VertexAttributes &v2, const VertexAttributes &v3, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
	{
		Eigen::Matrix<double, 3, 4> p;
		p.row(0) = v1.position.array() / v1.position[3];
		p.row(1) = v2.position.array() / v2.position[3];
		p.row(2) = v3.position.array() / v3.position[3];

		p.col(0) = ((p.col(0).array() + 1.0) / 2.0) * frameBuffer.rows();
		p.col(1) = ((p.col(1).array() + 1.0) / 2.0) * frameBuffer.cols();

		int lx = std::floor(p.col(0).minCoeff());
		int ly = std::floor(p.col(1).minCoeff());
		int ux = std::ceil(p.col(0).maxCoeff());
		int uy = std::ceil(p.col(1).maxCoeff());

		lx = std::min(std::max(lx, int(0)), int(frameBuffer.rows() - 1));
		ly = std::min(std::max(ly, int(0)), int(frameBuffer.cols() - 1));
		ux = std::max(std::min(ux, int(frameBuffer.rows() - 1)), int(0));
		uy = std::max(std::min(uy, int(frameBuffer.cols() - 1)), int(0));

		Eigen::Matrix3d A;
		A.col(0) = p.row(0).segment(0, 3);
		A.col(1) = p.row(1).segment(0, 3);
		A.col(2) = p.row(2).segment(0, 3);
		A.row(2) << 1.0, 1.0, 1.0;

		Eigen::Matrix3d Ai = A.inverse();

		for (unsigned i = lx; i <= ux; i++)
		{
			for (unsigned j = ly; j <= uy; j++)
			{

				Eigen::Vector3d pixel(i + 0.5, j + 0.5, 1);
				Eigen::Vector3d b = Ai * pixel;
				if (b.minCoeff() >= 0)
				{
					VertexAttributes va = VertexAttributes::interpolate(v1, v2, v3, b[0], b[1], b[2]);

					if (va.position[2] >= -1 && va.position[2] <= 1)
					{
						FragmentAttributes frag = program.FragmentShader(va, uniform);
						frameBuffer(i, j) = program.BlendingShader(frag, frameBuffer(i, j));
					}
				}
			}
		}
	}

	void rasterize_triangles(const Program &program, const UniformAttributes &uniform, const std::vector<VertexAttributes> &vertices, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
	{
		std::vector<VertexAttributes> v(vertices.size());
		for (unsigned i = 0; i < vertices.size(); i++)
			v[i] = program.VertexShader(vertices[i], uniform);

		for (unsigned i = 0; i < vertices.size() / 3; i++)
			rasterize_triangle(program, uniform, v[i * 3 + 0], v[i * 3 + 1], v[i * 3 + 2], frameBuffer);
	}

	void rasterize_line(const Program &program, const UniformAttributes &uniform, const VertexAttributes &v1, const VertexAttributes &v2, double line_thickness, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
	{
		Eigen::Matrix<double, 2, 4> p;
		p.row(0) = v1.position.array() / v1.position[3];
		p.row(1) = v2.position.array() / v2.position[3];

		p.col(0) = ((p.col(0).array() + 1.0) / 2.0) * frameBuffer.rows();
		p.col(1) = ((p.col(1).array() + 1.0) / 2.0) * frameBuffer.cols();

		int lx = std::floor(p.col(0).minCoeff() - line_thickness);
		int ly = std::floor(p.col(1).minCoeff() - line_thickness);
		int ux = std::ceil(p.col(0).maxCoeff() + line_thickness);
		int uy = std::ceil(p.col(1).maxCoeff() + line_thickness);

		lx = std::min(std::max(lx, int(0)), int(frameBuffer.rows() - 1));
		ly = std::min(std::max(ly, int(0)), int(frameBuffer.cols() - 1));
		ux = std::max(std::min(ux, int(frameBuffer.rows() - 1)), int(0));
		uy = std::max(std::min(uy, int(frameBuffer.cols() - 1)), int(0));

		Eigen::Vector2f l1(p(0, 0), p(0, 1));
		Eigen::Vector2f l2(p(1, 0), p(1, 1));

		double t = -1;
		double ll = (l1 - l2).squaredNorm();

		for (unsigned i = lx; i <= ux; i++)
		{
			for (unsigned j = ly; j <= uy; j++)
			{

				Eigen::Vector2f pixel(i + 0.5, j + 0.5);

				if (ll == 0.0)
					t = 0;
				else
				{
					t = (pixel - l1).dot(l2 - l1) / ll;
					t = std::fmax(0, std::fmin(1, t));
				}

				Eigen::Vector2f pixel_p = l1 + t * (l2 - l1);

				if ((pixel - pixel_p).squaredNorm() < (line_thickness * line_thickness))
				{
					VertexAttributes va = VertexAttributes::interpolate(v1, v2, v1, 1 - t, t, 0);
					FragmentAttributes frag = program.FragmentShader(va, uniform);
					frameBuffer(i, j) = program.BlendingShader(frag, frameBuffer(i, j));
				}
			}
		}
	}

	void rasterize_lines(const Program &program, const UniformAttributes &uniform, const std::vector<VertexAttributes> &vertices, double line_thickness, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
	{
		std::vector<VertexAttributes> v(vertices.size());
		for (unsigned i = 0; i < vertices.size(); i++)
			v[i] = program.VertexShader(vertices[i], uniform);

		for (unsigned i = 0; i < vertices.size() / 2; i++)
			rasterize_line(program, uniform, v[i * 2 + 0], v[i * 2 + 1], line_thickness, frameBuffer);
	}

	void framebuffer_to_uint8(const Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer, std::vector<uint8_t> &image)
	{
		const int w = frameBuffer.rows();
		const int h = frameBuffer.cols();
		const int comp = 4;
		const int stride_in_bytes = w * comp;
		image.resize(w * h * comp, 0);

		for (unsigned wi = 0; wi < w; ++wi)
		{
			for (unsigned hi = 0; hi < h; ++hi)
			{
				unsigned hif = h - 1 - hi;
				image[(hi * w * 4) + (wi * 4) + 0] = frameBuffer(wi, hif).color[0];
				image[(hi * w * 4) + (wi * 4) + 1] = frameBuffer(wi, hif).color[1];
				image[(hi * w * 4) + (wi * 4) + 2] = frameBuffer(wi, hif).color[2];
				image[(hi * w * 4) + (wi * 4) + 3] = frameBuffer(wi, hif).color[3];
			}
		}
	}
}