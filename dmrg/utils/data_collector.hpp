
#ifndef DATA_COLLECTOR_HPP_
#define DATA_COLLECTOR_HPP_

#include <string>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <map>

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#endif

#ifdef ENABLE_DATACOLLECTORS

#define DCOLLECTOR_CREATE(var, type, name) DataCollector<type> var(name);
#define DCOLLECTOR_GROUP(var, keyname) var.set_key(keyname);
#define DCOLLECTOR_ADD(var, value) var.add_data(value);
#define DCOLLECTOR_ADD_AT(var, keyname, value) var.add_data(keyname, value);
#define DCOLLECTOR_SAVE(var, ar, path) ar << alps::make_pvp(path, var);
#define DCOLLECTOR_SAVE_TO_FILE(var, fname, path)                           \
{                                                                           \
    alps::hdf5::archive h5ar_dcollector(fname, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE); \
    h5ar_dcollector << alps::make_pvp(path, var);                           \
}

#else

#define DCOLLECTOR_CREATE(var, type, name)
#define DCOLLECTOR_GROUP(var, keyname)
#define DCOLLECTOR_ADD(var, value)
#define DCOLLECTOR_ADD_AT(var, keyname, value)
#define DCOLLECTOR_SAVE(var, ar, path)
#define DCOLLECTOR_SAVE_TO_FILE(var, fname, path) 

#endif

template <class T>
class DataCollector
{
public:

	DataCollector(std::string const & name) : name_(name), active_key("none") {}

	std::string name() const {return name_;}

	void set_key (std::string const & key)
	{
		active_key = key;
	}

	void add_data (const T& val)
	{
		data[active_key].push_back(val);
	}
	void add_data (std::string const & key, const T& val)
	{
		data[key].push_back(val);
	}

#ifdef HAVE_ALPS_HDF5
	void save(alps::hdf5::archive & ar) const
	{
		if (data.size() == 1) {
			ar << alps::make_pvp("mean/value", data.begin()->second);
		} else if (data.size() > 1) {
			std::vector<std::string> keys;
            std::vector<std::vector<T> > values;
			for (typename std::map<std::string, std::vector<T> >::const_iterator it = data.begin();
				it != data.end();
				it++)
			{
                keys.push_back(it->first);
                values.push_back(it->second);
			}
			ar << alps::make_pvp("mean/value", values);
			ar << alps::make_pvp("labels", keys);
		}
	}
#endif

private:
	std::map<std::string, std::vector<T> > data;
	std::string active_key;
	std::string name_;
};


#endif /* DATA_COLLECTOR_HPP_ */
