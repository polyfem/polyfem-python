#pragma once

#include <polyfem/io/ResourceIO.hpp>
#include <polyfem/mesh/MeshLoader.hpp>
#include <stdexcept>

// Expose one array mesh through the backend's existing resource interface.
class ArrayMeshResources final : public polyfem::io::ResourceIO
{
public:
  ArrayMeshResources(const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &cells)
      : data_(std::make_shared<Data>(Data{vertices, cells})),
        files_(std::make_unique<polyfem::io::FileSystemIO>(std::filesystem::current_path())) {}

  static std::string mesh_path() { return "/__polyfem_python_array_mesh__"; }

  std::unique_ptr<const ResourceIO> with_root(const std::string &root) const override
  {
    return std::unique_ptr<const ResourceIO>(new ArrayMeshResources(data_, files_->with_root(root)));
  }
  std::string canonical_path(const std::string &path) const override
  { return owned(path) ? path : files_->canonical_path(path); }
  bool exists(const std::string &path) const override
  {
    if (!owned(path)) return files_->exists(path);
    return path == mesh_path() || path == mesh_path() + "/vertices" || path == mesh_path() + "/cells";
  }
  bool is_group(const std::string &path) const override
  { return owned(path) ? path == mesh_path() : files_->is_group(path); }
  std::vector<std::string> list(const std::string &path) const override
  { return owned(path) ? std::vector<std::string>{mesh_path() + "/vertices", mesh_path() + "/cells"} : files_->list(path); }
  Eigen::MatrixXd read_matrix(const std::string &path) const override
  { if (path == mesh_path() + "/vertices") return data_->vertices; return files_->read_matrix(path); }
  Eigen::MatrixXi read_int_matrix(const std::string &path) const override
  { if (path == mesh_path() + "/cells") return data_->cells; return files_->read_int_matrix(path); }
  bool has_attribute(const std::string &path, const std::string &name) const override
  {
    if (!owned(path)) return files_->has_attribute(path, name);
    return path == mesh_path() && (name == "schema_version" || name == "mesh_type" || name == "dimension");
  }
  long read_integer_attribute(const std::string &path, const std::string &name) const override
  {
    if (!owned(path)) return files_->read_integer_attribute(path, name);
    if (name == "schema_version") return polyfem::mesh::MESH_SCHEMA_VERSION;
    if (name == "dimension") return data_->vertices.cols();
    throw std::runtime_error("Unknown array mesh attribute: " + name);
  }
  std::string read_string_attribute(const std::string &path, const std::string &name) const override
  {
    if (!owned(path)) return files_->read_string_attribute(path, name);
    if (name == "mesh_type") return "fem";
    throw std::runtime_error("Unknown array mesh attribute: " + name);
  }
  const std::filesystem::path &host_directory() const override { return files_->host_directory(); }
  std::unique_ptr<std::istream> open(const std::string &path, bool binary) const override { return files_->open(path, binary); }
  std::vector<double> read_double_vector(const std::string &path) const override { return files_->read_double_vector(path); }
  std::vector<int> read_int_vector(const std::string &path) const override { return files_->read_int_vector(path); }
  std::vector<long> read_long_vector(const std::string &path) const override { return files_->read_long_vector(path); }
  std::array<long, 2> read_shape_attribute(const std::string &path, const std::string &name) const override { return files_->read_shape_attribute(path, name); }
  std::filesystem::path materialize(const std::string &path) const override { return files_->materialize(path); }
  std::string describe(const std::string &path) const override { return files_->describe(path); }
private:
  struct Data { Eigen::MatrixXd vertices; Eigen::MatrixXi cells; };
  std::shared_ptr<const Data> data_;
  std::unique_ptr<const ResourceIO> files_;
  static bool owned(const std::string &path)
  { return path == mesh_path() || path.rfind(mesh_path() + "/", 0) == 0; }
  ArrayMeshResources(std::shared_ptr<const Data> data, std::unique_ptr<const ResourceIO> files)
      : data_(std::move(data)), files_(std::move(files)) {}
};
