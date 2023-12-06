import './App.css';
import { CCP4i2MoorhenContainer } from './wrapper/CCP4i2MoorhenContainer';
import { MoorhenReduxProvider } from 'moorhen';

import {
  createBrowserRouter, RouterProvider, useParams,
} from "react-router-dom";

function App() {
  return (<RouterProvider router={router} />
  );
}

const router = createBrowserRouter([
  {
    path: "/CCP4i2Moorhen/",
    element: (
      <CCP4i2Moorhen />
    ),
  },
  {
    path: "/",
    element: (
      <CCP4i2Moorhen />
    ),
  },
]);


function CCP4i2Moorhen() {
  const mySearchParams = new URLSearchParams(window.location.search.slice(1));
  //Provide a default for the port number...DANGER only relevant for development
  const searchDict = {CCP4i2Port:43434}

  for (const [key, value] of mySearchParams.entries()) {
    searchDict[key] = decodeURIComponent(value)
  }
  return <div className="App">
    <MoorhenReduxProvider>
      <CCP4i2MoorhenContainer {...searchDict} urlRoot={`http://127.0.0.1:${searchDict.CCP4i2Port}/database`}/>
    </MoorhenReduxProvider>
  </div>
}
export default App;
