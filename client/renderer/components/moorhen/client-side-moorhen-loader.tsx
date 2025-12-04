import dynamic from "next/dynamic";
import { MoorhenWrapperProps } from "./moorhen-wrapper";

const MoorhenClient = dynamic(() => import("./client-side-moorhen-component"), {
  ssr: false,
});

export default function DynamicMoorhenClient(props: MoorhenWrapperProps) {
  return <MoorhenClient {...props} />;
}
